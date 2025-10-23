import json
import logging
from dataclasses import asdict
from pathlib import Path
from typing import Optional, Sequence

from skbio.alignment import StripedSmithWaterman  # type: ignore

from riot_na.alignment.alignment_utils import infer_reading_frame, translate
from riot_na.alignment.gene_aligner import (
    GeneAlignerAA,
    create_aa_c_gene_aligner,
    create_aa_j_gene_aligner,
    create_aa_v_gene_aligner,
)
from riot_na.alignment.segment_gene_aligner import (
    AA_ALIGNER_PARAMS,
    SegmentGeneAlignerAA,
)
from riot_na.alignment.segment_gene_aligner import (
    create_aa_v_gene_aligner as create_segment_aa_v_gene_aligner,
)
from riot_na.alignment.skbio_alignment import align_aa
from riot_na.config import GENE_DB_DIR
from riot_na.data.model import (
    AlignmentEntryAA,
    AlignmentsAA,
    AlignmentSegment,
    AlignmentsNT,
    InternalAlignmentEntryAA,
    Locus,
    Organism,
    TranslatedAlignmentsNT,
)

logger = logging.getLogger(__name__)


class VJCAlignerAA:
    def __init__(
        self,
        v_aligner: SegmentGeneAlignerAA | GeneAlignerAA,
        j_aligners: dict[Organism, dict[Locus, GeneAlignerAA]],
        c_aligners: Optional[dict[Organism, dict[Locus, GeneAlignerAA]]] = None,
    ) -> None:
        self.v_aligner = v_aligner
        self.j_aligners = j_aligners
        self.c_aligners = c_aligners

    def produce_aa_alignments(self, query: str) -> list[AlignmentsAA]:
        results = []

        v_segment_alignments: Sequence[AlignmentSegment] = self.v_aligner.align(query)

        for v_aa_segment_alignment in v_segment_alignments:  # pylint: disable=too-many-nested-blocks
            try:
                result = self._produce_vjc_aa_alignment(v_aa_segment_alignment)
                results.append(result)
            except ValueError as e:
                logger.error(
                    "Error aligning segment %s due to %s",
                    (
                        v_aa_segment_alignment.best_alignment.q_seq
                        if v_aa_segment_alignment.best_alignment
                        else "unknown sequence"
                    ),
                    e,
                )
        return results

    def _produce_vjc_aa_alignment(self, v_aa_segment_alignment: AlignmentSegment) -> AlignmentsAA:
        result = AlignmentsAA()
        result.segment_start = v_aa_segment_alignment.segment_start
        result.segment_end = v_aa_segment_alignment.segment_end

        v_aa_alignment = v_aa_segment_alignment.best_alignment
        if v_aa_alignment is None:
            return result

        result.v = v_aa_alignment

        query_j = v_aa_alignment.q_seq[v_aa_alignment.q_end :]
        if query_j:
            # align j genes
            species = v_aa_alignment.species
            locus = v_aa_alignment.locus

            j_aligner = self.j_aligners[species][locus]
            j_aa_alignment = j_aligner.align(query_j)[0].best_alignment

            if j_aa_alignment is None:
                return result

            result.j = j_aa_alignment

            if self.c_aligners:
                query_c = query_j[j_aa_alignment.q_end :]
                if query_c:
                    c_aligner = self.c_aligners[species][locus]
                    assert c_aligner is not None
                    c_aa_alignment = c_aligner.align(query_c)[0].best_alignment
                    if c_aa_alignment is not None:
                        result.c = c_aa_alignment

        return result


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
            # Note: SegmentGeneAlignerAA doesn't have gene_lookup attribute like GeneAlignerAA
            # We'll need to access the gene sequences directly from the aligner
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
                alignments, aligned_sequence_aa, species, locus, extended_j_aa_alignment, v_aa_alignment.q_end
            )

            if self.c_aligners is not None:
                c_aln_length_threshold = self.c_aligners[species][locus].alignment_length_threshold or 0
                if c_aa_alignment is not None and (  # pylint: disable=too-many-boolean-expressions
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
                    species=v_alignment.species,
                    locus=v_alignment.locus,
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
        locus: Locus,
        extended_j_aa_alignment: AlignmentEntryAA | None,
        j_offset: int,
    ) -> tuple[str, str | None, str | None, InternalAlignmentEntryAA | None]:
        if alignments.c and extended_j_aa_alignment:
            masked_aa_for_c = aligned_sequence_aa[extended_j_aa_alignment.q_end + j_offset :]

            c_gene_id = alignments.c.target_id

            assert self.c_aligners
            c_gene = self.c_aligners[species][locus].gene_lookup[alignments.c.lookup_gene_id()]
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


def create_vjc_aligner_aa(
    allowed_species: Optional[list[Organism]] = None, db_dir: Path = GENE_DB_DIR, use_segment_aligner: bool = False
) -> VJCAlignerAA:
    if not allowed_species:
        allowed_species = [Organism.HOMO_SAPIENS, Organism.MUS_MUSCULUS, Organism.VICUGNA_PACOS]

    aa_genes_dir = db_dir / "gene_db" / "aa_genes_deduplicated"
    v_aligner = (
        create_segment_aa_v_gene_aligner(allowed_species=allowed_species, aa_genes_dir=aa_genes_dir)
        if use_segment_aligner
        else create_aa_v_gene_aligner(allowed_species=allowed_species, aa_genes_dir=aa_genes_dir)
    )

    j_aligners: dict[Organism, dict[Locus, GeneAlignerAA]] = {}

    for organism in allowed_species:
        for locus in Locus:
            if organism == Organism.VICUGNA_PACOS and locus != Locus.IGH:
                continue

            j_aligner = create_aa_j_gene_aligner(organism=organism, locus=locus, aa_genes_dir=aa_genes_dir)
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

            c_aligner = create_aa_c_gene_aligner(organism=c_organism, locus=locus, aa_genes_dir=aa_genes_dir)
            c_organism_aligner[locus] = c_aligner
        c_aligners[organism] = c_organism_aligner

    vjc_aligner = VJCAlignerAA(v_aligner, j_aligners, c_aligners)

    return vjc_aligner


def create_vjc_alignment_translator_aa(
    allowed_species: Optional[list[Organism]] = None, db_dir: Path = GENE_DB_DIR
) -> VJCAlignmentTranslatorAA:
    if not allowed_species:
        allowed_species = [Organism.HOMO_SAPIENS, Organism.MUS_MUSCULUS, Organism.VICUGNA_PACOS]

    aa_genes_dir = db_dir / "gene_db" / "aa_genes"
    v_aligner = create_aa_v_gene_aligner(allowed_species=allowed_species, aa_genes_dir=aa_genes_dir)

    j_aligners: dict[Organism, dict[Locus, GeneAlignerAA]] = {}

    for organism in allowed_species:
        for locus in Locus:
            if organism == Organism.VICUGNA_PACOS and locus != Locus.IGH:
                continue

            j_aligner = create_aa_j_gene_aligner(organism=organism, locus=locus, aa_genes_dir=aa_genes_dir)
            organism_aligners = j_aligners.get(organism, {})
            organism_aligners[locus] = j_aligner
            j_aligners[organism] = organism_aligners

    c_aligners: dict[Organism, dict[Locus, GeneAlignerAA]] = {}

    for organism in allowed_species:
        c_organism_aligner = c_aligners.get(organism, {})

        if organism not in [Organism.HOMO_SAPIENS, Organism.VICUGNA_PACOS, Organism.CUSTOM]:
            c_organism = Organism.HOMO_SAPIENS
        else:
            c_organism = organism

        for locus in Locus:

            # at the moment we only have heavy genes
            if organism == Organism.VICUGNA_PACOS and locus != Locus.IGH:
                continue

            c_aligner = create_aa_c_gene_aligner(organism=c_organism, locus=locus, aa_genes_dir=aa_genes_dir)
            c_organism_aligner[locus] = c_aligner
        c_aligners[organism] = c_organism_aligner

    vj_aligner = VJCAlignmentTranslatorAA(v_aligner, j_aligners, c_aligners)

    return vj_aligner


if __name__ == "__main__":
    # case 1 detect organism

    # given query
    SAMPLE_QUERY = "QVQLQQWGAGLLKPSETLSLTCAVFGGSFSGYYWSWIRQPPGKGLEWIGEINHRGNTNDNPSLKSRVTISVDTSKNQFALKLSSVTAADTAVYYCARERGYTYGNFDHWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"

    vjc_alnr = create_vjc_aligner_aa(allowed_species=[Organism.HOMO_SAPIENS])

    local_aa_alignments = vjc_alnr.produce_aa_alignments(SAMPLE_QUERY)
    print(json.dumps([asdict(aa_alignment) for aa_alignment in local_aa_alignments], indent=4))
