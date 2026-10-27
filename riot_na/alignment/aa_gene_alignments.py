import json
import logging
from dataclasses import asdict
from pathlib import Path
from typing import Optional

from ssw_aligner import StripedSmithWaterman

from riot_na.alignment.alignment_utils import align_c, infer_reading_frame
from riot_na.alignment.gene_aligner import (
    AA_ALIGNER_PARAMS,
    GENE_SELECTION_SPECIES_ORDER,
    GeneAlignerAA,
    create_aa_c_gene_aligner,
    create_aa_j_gene_aligner,
    create_aa_v_gene_aligner,
)
from riot_na.alignment.skbio_alignment import align_aa
from riot_na.alignment.translate import translate
from riot_na.config import GENE_DB_DIR
from riot_na.data.model import (
    AlignmentEntryAA,
    AlignmentsAA,
    AlignmentsNT,
    InternalAlignmentEntryAA,
    Locus,
    Organism,
    TranslatedAlignmentsNT,
)

logger = logging.getLogger(__name__)


def align_jc_aa(
    v_aa_alignment: AlignmentEntryAA,
    j_aligners: dict[Organism, dict[Locus, GeneAlignerAA]],
    c_aligners: Optional[dict[Organism, dict[Locus, GeneAlignerAA]]],
) -> AlignmentsAA:
    """Align J and C genes to the part of the query following the V alignment."""
    result = AlignmentsAA(v=v_aa_alignment)

    query_j = v_aa_alignment.q_seq[v_aa_alignment.q_end :]
    if query_j:
        # align j genes
        species = v_aa_alignment.species
        locus = v_aa_alignment.locus

        j_aligner = j_aligners[species][locus]
        j_aa_alignment = j_aligner.align(query_j)

        # C follows J or V, if J was not found
        query_c = query_j
        if j_aa_alignment is not None:
            result.j = j_aa_alignment
            query_c = query_j[j_aa_alignment.q_end :]

        if c_aligners and query_c:
            result.c = align_c(c_aligners[species], locus, query_c, lambda c_aligner: c_aligner.align(query_c))

    return result


class VJCAlignerAA:
    def __init__(
        self,
        v_aligner: GeneAlignerAA,
        j_aligners: dict[Organism, dict[Locus, GeneAlignerAA]],
        c_aligners: Optional[dict[Organism, dict[Locus, GeneAlignerAA]]] = None,
    ) -> None:
        self.v_aligner = v_aligner
        self.j_aligners = j_aligners
        self.c_aligners = c_aligners

    def produce_aa_alignments(self, query: str) -> AlignmentsAA:
        v_aa_alignment = self.v_aligner.align(query)
        if v_aa_alignment is None:
            return AlignmentsAA()

        try:
            return align_jc_aa(v_aa_alignment, self.j_aligners, self.c_aligners)
        except ValueError as e:
            logger.error("Error aligning %s due to %s", v_aa_alignment.q_seq, e)
            return AlignmentsAA()

    def produce_aa_alignments_for_genes(self, query: str, v_gene: str, j_gene: str) -> AlignmentsAA:
        """
        Align query to the given V and J germlines, skipping prefiltering.

        Germlines are given by their exact names, e.g. IGHV3-7*01 and IGHJ4*02, and looked up in the species
        the aligner was created for, human first. The J germline is aligned to the part of the query left over
        by the V alignment, so the J locus follows from the V germline. C genes are not assigned.
        """
        result = AlignmentsAA()

        v_aa_alignments = self.v_aligner.align_to_genes(query, [v_gene])
        if not v_aa_alignments:
            return result

        v_aa_alignment = v_aa_alignments[0]
        result.v = v_aa_alignment

        query_j = v_aa_alignment.q_seq[v_aa_alignment.q_end :]
        if query_j:
            j_aligner = self._get_j_aligner(v_aa_alignment.locus, j_gene)
            j_aa_alignments = j_aligner.align_to_genes(query_j, [j_gene])
            if j_aa_alignments:
                result.j = j_aa_alignments[0]

        return result

    def assert_genes_exist(self, v_gene: str, j_gene: str) -> None:
        """Fail fast on unknown germline names, without aligning anything."""
        v_germline = self.v_aligner.get_gene(v_gene)
        self._get_j_aligner(v_germline.locus, j_gene)

    def _get_j_aligner(self, locus: Locus, j_gene: str) -> GeneAlignerAA:
        """J aligner of the locus with the given germline, human first."""
        for species in GENE_SELECTION_SPECIES_ORDER:
            j_aligner = self.j_aligners.get(species, {}).get(locus)
            if j_aligner is not None and j_aligner.find_gene(j_gene) is not None:
                return j_aligner
        raise AssertionError(f"Gene {j_gene} not found in {locus.value} J germline database")


class VJCAlignmentTranslatorAA:
    def __init__(
        self,
        v_aligner: GeneAlignerAA,
        j_aligners: dict[Organism, dict[Locus, GeneAlignerAA]],
        c_aligners: Optional[dict[Organism, dict[Locus, GeneAlignerAA]]] = None,
    ) -> None:
        self.v_aligner = v_aligner
        self.j_aligners = j_aligners
        self.c_aligners = c_aligners

    def translate_nt_alignments(  # pylint: disable=too-many-statements
        self, alignments: AlignmentsNT, extend_alignment: bool = True
    ) -> Optional[TranslatedAlignmentsNT]:
        if alignments.v:
            # Query might have been reverse complemented during V alignment.
            query_sequence = alignments.v.q_seq
            v_alignment = alignments.v
            assert v_alignment.reading_frame is not None
            reading_frame = infer_reading_frame(t_start=v_alignment.t_start, t_frame=v_alignment.reading_frame)

            # translate sequence to aa
            alignment_start = alignments.v.q_start
            alignment_end = alignments.v.q_end

            if alignments.c:
                alignment_end = alignments.c.q_end
            elif alignments.j:
                alignment_end = alignments.j.q_end

            aligned_sequence = query_sequence[alignment_start:alignment_end]
            aligned_sequence_aa = translate(aligned_sequence, reading_frame)

            aligned_sequence_no_c = query_sequence[
                alignment_start : alignments.j.q_end if alignments.j else alignments.v.q_end
            ]
            alignes_sequence_no_c_aa = translate(aligned_sequence_no_c, reading_frame)

            species = v_alignment.species
            locus = v_alignment.locus
            v_gene = self.v_aligner.gene_lookup.get(v_alignment.lookup_gene_id())
            v_gene_aa = v_gene.sequence if v_gene else ""

            locus = v_alignment.locus

            j_gene_id = alignments.j.target_id if alignments.j else None
            j_gene_aa = None
            if j_gene_id:
                assert alignments.j
                oragnism_locus_j_aligner = self.j_aligners[species][locus]
                assert oragnism_locus_j_aligner is not None
                j_gene = oragnism_locus_j_aligner.gene_lookup.get(alignments.j.lookup_gene_id())
                j_gene_aa = j_gene.sequence if j_gene else ""

            v_aligner = StripedSmithWaterman(aligned_sequence_aa, **AA_ALIGNER_PARAMS)
            v_aa_alignment = align_aa(
                v_aligner,
                aligned_sequence_aa,
                alignments.v.target_id,
                v_gene_aa,
                calculate_score=False,
                extend_n_terminus=extend_alignment,
            )
            v_aln_length_threshold = self.v_aligner.alignment_length_threshold or 0
            assert (
                v_aa_alignment.t_start != -1
                and v_aa_alignment.t_end != -1
                and v_aa_alignment.q_start != -1
                and v_aa_alignment.q_end != -1
                and v_aa_alignment.q_end - v_aa_alignment.q_start >= v_aln_length_threshold
            ), "Query translation failed: could not align sequence_aa to germline_aa"

            masked_aa = alignes_sequence_no_c_aa[v_aa_alignment.q_end :]
            j_aligner = StripedSmithWaterman(masked_aa, **AA_ALIGNER_PARAMS)
            j_aa_alignment = (
                align_aa(
                    j_aligner,
                    masked_aa,
                    j_gene_id,
                    j_gene_aa,
                    calculate_score=False,
                    extend_c_terminus=extend_alignment,
                )
                if j_gene_id and j_gene_aa
                else None
            )

            j_aln_length_threshold = self.j_aligners[species][locus].alignment_length_threshold or 0
            if j_aa_alignment is not None and (  # pylint: disable=too-many-boolean-expressions
                j_aa_alignment.q_start == -1
                or j_aa_alignment.q_end == -1
                or j_aa_alignment.t_start == -1
                or j_aa_alignment.t_end == -1
                or j_aa_alignment.q_end - j_aa_alignment.q_start < j_aln_length_threshold
            ):
                j_aa_alignment = None

            extended_v_aa_alignment = AlignmentEntryAA(
                target_id=v_alignment.target_id,
                alignment_score=v_aa_alignment.alignment_score,
                seq_identity=v_aa_alignment.seq_identity,
                e_value=v_aa_alignment.e_value,
                q_start=v_aa_alignment.q_start,
                q_end=v_aa_alignment.q_end,
                t_start=v_aa_alignment.t_start,
                t_end=v_aa_alignment.t_end,
                cigar=v_aa_alignment.cigar,
                species=v_alignment.species,
                locus=v_alignment.locus,
                q_seq=aligned_sequence_aa,
                t_seq=v_gene_aa,
            )

            if j_aa_alignment:
                assert j_gene_aa
                assert j_gene_id

                extended_j_aa_alignment = AlignmentEntryAA(
                    target_id=j_gene_id,
                    alignment_score=j_aa_alignment.alignment_score,
                    seq_identity=j_aa_alignment.seq_identity,
                    e_value=j_aa_alignment.e_value,
                    q_start=j_aa_alignment.q_start,
                    q_end=j_aa_alignment.q_end,
                    t_start=j_aa_alignment.t_start,
                    t_end=j_aa_alignment.t_end,
                    cigar=j_aa_alignment.cigar,
                    species=v_alignment.species,
                    locus=v_alignment.locus,
                    q_seq=masked_aa,
                    t_seq=j_gene_aa,
                )
            else:
                extended_j_aa_alignment = None

            masked_aa_for_c, c_gene_id, c_gene_aa, c_aa_alignment = self.generate_c_aa_alignment(
                alignments, aligned_sequence_aa, species, extended_j_aa_alignment, v_aa_alignment.q_end
            )

            if self.c_aligners is not None and c_aa_alignment is not None:
                assert alignments.c
                c_aln_length_threshold = self.c_aligners[species][alignments.c.locus].alignment_length_threshold or 0
                if (  # pylint: disable=too-many-boolean-expressions
                    c_aa_alignment.q_start == -1
                    or c_aa_alignment.q_end == -1
                    or c_aa_alignment.t_start == -1
                    or c_aa_alignment.t_end == -1
                    or c_aa_alignment.q_end - c_aa_alignment.q_start < c_aln_length_threshold
                ):
                    c_aa_alignment = None

            if c_aa_alignment:
                assert c_gene_aa
                assert c_gene_id
                assert alignments.c

                # C gene can be of another locus and species than V gene
                extended_c_aa_alignment = AlignmentEntryAA(
                    target_id=c_gene_id,
                    alignment_score=c_aa_alignment.alignment_score,
                    seq_identity=c_aa_alignment.seq_identity,
                    e_value=c_aa_alignment.e_value,
                    q_start=c_aa_alignment.q_start,
                    q_end=c_aa_alignment.q_end,
                    t_start=c_aa_alignment.t_start,
                    t_end=c_aa_alignment.t_end,
                    cigar=c_aa_alignment.cigar,
                    species=alignments.c.species,
                    locus=alignments.c.locus,
                    q_seq=masked_aa_for_c,
                    t_seq=c_gene_aa,
                )
            else:
                extended_c_aa_alignment = None

            aa_alignments = AlignmentsAA(
                v=extended_v_aa_alignment, j=extended_j_aa_alignment, c=extended_c_aa_alignment
            )

            return TranslatedAlignmentsNT(
                translated_query=aligned_sequence_aa, reading_frame=reading_frame, aa_alignments=aa_alignments
            )

        return None

    def generate_j_aa_alignment(
        self,
        aligned_sequence_aa: str,
        j_gene_id: Optional[str],
        j_gene_aa: Optional[str],
        v_aa_alignment: InternalAlignmentEntryAA,
    ) -> tuple[str, InternalAlignmentEntryAA | None]:
        masked_aa = aligned_sequence_aa[v_aa_alignment.q_end :]
        j_aligner = StripedSmithWaterman(masked_aa, **AA_ALIGNER_PARAMS)
        j_aa_alignment = (
            align_aa(j_aligner, masked_aa, j_gene_id, j_gene_aa, calculate_score=False)
            if j_gene_id and j_gene_aa
            else None
        )

        if j_aa_alignment is not None and (
            j_aa_alignment.q_start == -1
            or j_aa_alignment.q_end == -1
            or j_aa_alignment.t_start == -1
            or j_aa_alignment.t_end == -1
        ):
            j_aa_alignment = None
        return masked_aa, j_aa_alignment

    def generate_c_aa_alignment(
        self,
        alignments: AlignmentsNT,
        aligned_sequence_aa: str,
        species: Organism,
        extended_j_aa_alignment: AlignmentEntryAA | None,
        j_offset: int,
    ) -> tuple[str, str | None, str | None, InternalAlignmentEntryAA | None]:
        if alignments.c:
            # C follows J or V, if J was not found
            c_offset = j_offset + extended_j_aa_alignment.q_end if extended_j_aa_alignment else j_offset
            masked_aa_for_c = aligned_sequence_aa[c_offset:]

            c_gene_id = alignments.c.target_id

            assert self.c_aligners
            c_gene = self.c_aligners[species][alignments.c.locus].gene_lookup[alignments.c.lookup_gene_id()]
            c_gene_aa = c_gene.sequence

            c_aligner = StripedSmithWaterman(masked_aa_for_c, **AA_ALIGNER_PARAMS)
            c_aa_alignment = align_aa(
                c_aligner, masked_aa_for_c, alignments.c.target_id, c_gene_aa, calculate_score=False
            )
        else:
            masked_aa_for_c = ""
            c_gene_id = None
            c_gene_aa = None
            c_aa_alignment = None
        return masked_aa_for_c, c_gene_id, c_gene_aa, c_aa_alignment


DEFAULT_SPECIES = (Organism.HOMO_SAPIENS, Organism.MUS_MUSCULUS, Organism.VICUGNA_PACOS)


def aa_genes_dir(db_dir: Path) -> Path:
    return db_dir / "gene_db" / "aa_genes_deduplicated"


def create_jc_aligners_aa(
    allowed_species: tuple[Organism, ...], genes_dir: Path
) -> tuple[dict[Organism, dict[Locus, GeneAlignerAA]], dict[Organism, dict[Locus, GeneAlignerAA]]]:
    """J and C aligners of the allowed species by locus."""
    j_aligners: dict[Organism, dict[Locus, GeneAlignerAA]] = {}

    for organism in allowed_species:
        for locus in Locus:
            if organism == Organism.VICUGNA_PACOS and locus != Locus.IGH:
                continue

            j_aligner = create_aa_j_gene_aligner(organism=organism, locus=locus, aa_genes_dir=genes_dir)
            organism_aligners = j_aligners.get(organism, {})
            organism_aligners[locus] = j_aligner
            j_aligners[organism] = organism_aligners

    c_aligners: dict[Organism, dict[Locus, GeneAlignerAA]] = {}

    for organism in allowed_species:
        if organism not in [Organism.HOMO_SAPIENS, Organism.VICUGNA_PACOS, Organism.CUSTOM]:
            c_organism = Organism.HOMO_SAPIENS
        else:
            c_organism = organism

        c_organism_aligner = c_aligners.get(organism, {})

        for locus in Locus:

            # at the moment we only have heavy genes
            if organism == Organism.VICUGNA_PACOS and locus != Locus.IGH:
                continue

            c_aligner = create_aa_c_gene_aligner(organism=c_organism, locus=locus, aa_genes_dir=genes_dir)
            c_organism_aligner[locus] = c_aligner
        c_aligners[organism] = c_organism_aligner

    return j_aligners, c_aligners


def create_vjc_aligner_aa(
    allowed_species: Optional[tuple[Organism, ...]] = None,
    db_dir: Path = GENE_DB_DIR,
) -> VJCAlignerAA:
    allowed_species = allowed_species or DEFAULT_SPECIES
    genes_dir = aa_genes_dir(db_dir)

    v_aligner = create_aa_v_gene_aligner(allowed_species=allowed_species, aa_genes_dir=genes_dir)
    j_aligners, c_aligners = create_jc_aligners_aa(allowed_species, genes_dir)

    return VJCAlignerAA(v_aligner, j_aligners, c_aligners)


def create_vjc_alignment_translator_aa(
    allowed_species: Optional[tuple[Organism, ...]] = None, db_dir: Path = GENE_DB_DIR
) -> VJCAlignmentTranslatorAA:
    allowed_species = allowed_species or DEFAULT_SPECIES
    genes_dir = db_dir / "gene_db" / "aa_genes"

    v_aligner = create_aa_v_gene_aligner(allowed_species=allowed_species, aa_genes_dir=genes_dir)
    j_aligners, c_aligners = create_jc_aligners_aa(allowed_species, genes_dir)

    return VJCAlignmentTranslatorAA(v_aligner, j_aligners, c_aligners)


if __name__ == "__main__":
    # case 1 detect organism

    # given query
    SAMPLE_QUERY = "QVQLQQWGAGLLKPSETLSLTCAVFGGSFSGYYWSWIRQPPGKGLEWIGEINHRGNTNDNPSLKSRVTISVDTSKNQFALKLSSVTAADTAVYYCARERGYTYGNFDHWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"

    vjc_alnr = create_vjc_aligner_aa(allowed_species=(Organism.HOMO_SAPIENS,))

    local_aa_alignments = vjc_alnr.produce_aa_alignments(SAMPLE_QUERY)
    print(json.dumps(asdict(local_aa_alignments), indent=4))
