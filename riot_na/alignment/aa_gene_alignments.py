import json
from dataclasses import asdict
from pathlib import Path
from typing import Optional

from skbio.alignment import StripedSmithWaterman  # type: ignore

from riot_na.alignment.alignment_utils import infer_reading_frame, translate
from riot_na.alignment.gene_aligner import (
    AA_ALIGNER_PARAMS,
    GeneAlignerAA,
    create_aa_j_gene_aligner,
    create_aa_v_gene_aligner,
)
from riot_na.alignment.skbio_alignment import align_aa
from riot_na.config import GENE_DB_DIR
from riot_na.data.model import (
    AlignmentEntryAA,
    AlignmentsAA,
    AlignmentsNT,
    Locus,
    Organism,
    TranslatedAlignmentsNT,
)


class VJAlignerAA:
    def __init__(
        self,
        v_aligner: GeneAlignerAA,
        j_aligners: dict[Organism, dict[Locus, GeneAlignerAA]],
    ) -> None:
        self.v_aligner = v_aligner
        self.j_aligners = j_aligners

    def produce_aa_alignments(self, query: str) -> AlignmentsAA:
        result = AlignmentsAA()

        v_aa_alignment = self.v_aligner.align(query)

        if v_aa_alignment is None:
            return result

        result.v = v_aa_alignment

        # mask alignments
        query_j = query[v_aa_alignment.q_end :]

        # align j genes
        species = v_aa_alignment.species
        locus = v_aa_alignment.locus

        j_aligner = self.j_aligners[species][locus]
        j_aa_alignment = j_aligner.align(query_j)

        if j_aa_alignment is None:
            return result

        assert j_aa_alignment is not None
        result.j = j_aa_alignment

        return result


class VJAlignmentTranslatorAA:
    def __init__(
        self,
        v_aligner: GeneAlignerAA,
        j_aligners: dict[Organism, dict[Locus, GeneAlignerAA]],
    ) -> None:
        self.v_aligner = v_aligner
        self.j_aligners = j_aligners

    def translate_nt_alignments(self, alignments: AlignmentsNT) -> Optional[TranslatedAlignmentsNT]:
        if alignments.v:
            # Query might have been reverse complemented during V alignment.
            query_sequence = alignments.v.q_seq
            v_alignment = alignments.v
            assert v_alignment.reading_frame is not None
            reading_frame = infer_reading_frame(t_start=v_alignment.t_start, t_frame=v_alignment.reading_frame)

            # translate sequence to aa
            alignment_start = alignments.v.q_start
            alignment_end = alignments.j.q_end if alignments.j else alignments.v.q_end

            aligned_sequence = query_sequence[alignment_start:alignment_end]
            aligned_sequence_aa = translate(aligned_sequence, reading_frame)

            species = v_alignment.species
            v_gene = self.v_aligner.gene_lookup[species + "|" + v_alignment.target_id]

            locus = v_alignment.locus
            v_gene_aa = v_gene.sequence

            j_gene_id = alignments.j.target_id if alignments.j else None
            j_gene = self.j_aligners[species][locus].gene_lookup[species + "|" + j_gene_id] if j_gene_id else None

            j_gene_aa = j_gene.sequence if j_gene else None

            v_aligner = StripedSmithWaterman(aligned_sequence_aa, **AA_ALIGNER_PARAMS)
            v_aa_alignment = align_aa(
                v_aligner, aligned_sequence_aa, alignments.v.target_id, v_gene_aa, calculate_score=False
            )
            assert (
                v_aa_alignment.t_start != -1
                and v_aa_alignment.t_end != -1
                and v_aa_alignment.q_start != -1
                and v_aa_alignment.q_end != -1
            ), "Query translation failed: could not align sequence_aa to germline_aa"

            masked_aa = aligned_sequence_aa[v_aa_alignment.q_end :]
            j_aligner = StripedSmithWaterman(masked_aa, **AA_ALIGNER_PARAMS)
            j_aa_alignment = (
                align_aa(j_aligner, masked_aa, j_gene_id, j_gene_aa, calculate_score=False)
                if j_gene_id and j_gene_aa
                else None
            )

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

            aa_alignments = AlignmentsAA(v=extended_v_aa_alignment, j=extended_j_aa_alignment)

            return TranslatedAlignmentsNT(
                translated_query=aligned_sequence_aa, reading_frame=reading_frame, aa_alignments=aa_alignments
            )

        return None


def create_vj_aligner_aa(allowed_species: Optional[list[Organism]] = None, db_dir: Path = GENE_DB_DIR) -> VJAlignerAA:
    if not allowed_species:
        allowed_species = list(Organism)

    aa_genes_dir = db_dir / "gene_db" / "aa_genes_deduplicated"
    v_aligner = create_aa_v_gene_aligner(allowed_species=allowed_species, aa_genes_dir=aa_genes_dir)

    j_aligners: dict[Organism, dict[Locus, GeneAlignerAA]] = {}

    for organism in allowed_species:
        for locus in Locus:
            j_aligner = create_aa_j_gene_aligner(organism=organism, locus=locus, aa_genes_dir=aa_genes_dir)
            organism_aligners = j_aligners.get(organism, {})
            organism_aligners[locus] = j_aligner
            j_aligners[organism] = organism_aligners

    vj_aligner = VJAlignerAA(v_aligner, j_aligners)

    return vj_aligner


def create_vj_alignment_translator_aa(
    allowed_species: Optional[list[Organism]] = None, db_dir: Path = GENE_DB_DIR
) -> VJAlignmentTranslatorAA:
    if not allowed_species:
        allowed_species = list(Organism)

    aa_genes_dir = db_dir / "gene_db" / "aa_genes"
    v_aligner = create_aa_v_gene_aligner(allowed_species=allowed_species, aa_genes_dir=aa_genes_dir)

    j_aligners: dict[Organism, dict[Locus, GeneAlignerAA]] = {}

    for organism in allowed_species:
        for locus in Locus:
            j_aligner = create_aa_j_gene_aligner(organism=organism, locus=locus, aa_genes_dir=aa_genes_dir)
            organism_aligners = j_aligners.get(organism, {})
            organism_aligners[locus] = j_aligner
            j_aligners[organism] = organism_aligners

    vj_aligner = VJAlignmentTranslatorAA(v_aligner, j_aligners)

    return vj_aligner


if __name__ == "__main__":
    # case 1 detect organism

    # given query
    SAMPLE_QUERY = (
        "LRRPCPSPALSLMAPSVVTTGTGSGSPQGRDWSGLGIYTIVGAPTTTPPSRVESPYQTPKNHFSLKLRSVTAVDTAVYYCARWGHFDTSGYFVVDYWGQGTLVTVSS"
    )

    vdj_alnr = create_vj_aligner_aa()

    nt_alignments = vdj_alnr.produce_aa_alignments(SAMPLE_QUERY)
    print(json.dumps(asdict(nt_alignments), indent=4))
