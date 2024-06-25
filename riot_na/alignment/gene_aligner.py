import os
from functools import partial
from pathlib import Path
from typing import Callable, Optional, Sequence

import blosum  # type: ignore
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from skbio.alignment import StripedSmithWaterman  # type: ignore

from riot_na.alignment.skbio_alignment import (
    align_aa,  # type: ignore  # pylint: disable=import-error
)
from riot_na.alignment.skbio_alignment import align
from riot_na.common.multi_species_prefiltering import MultiSpeciesPrefiltering
from riot_na.config import GENE_DB_DIR
from riot_na.data.model import (
    AlignmentEntryAA,
    AlignmentEntryNT,
    Gene,
    GeneAA,
    GermlineGene,
    Locus,
    Organism,
    SpeciesPrefilteringResult,
)

ALIGNER_PARAMS = {"match_score": 1, "mismatch_score": -1, "gap_open_penalty": 4, "gap_extend_penalty": 1}


class GeneAligner:
    def __init__(
        self,
        genes: list[Gene],
        kmer_size: int,
        distance_threshold: int,
        top_n: int,
        modulo_n: int,
        e_value_threshold: Optional[float] = None,
        alignment_length_threshold: Optional[int] = None,
    ) -> None:
        self.prefiltering = MultiSpeciesPrefiltering(genes, kmer_size, distance_threshold, top_n, modulo_n)

        self.gene_lookup = dict(map(lambda gene: (gene.species + "|" + gene.name, gene), genes))
        self.db_length = sum(len(gene.sequence) for gene in genes)

        self.e_value_threshold = e_value_threshold
        self.alignment_length_threshold = alignment_length_threshold

    def _prefilter(self, query: str, both_strains: bool) -> SpeciesPrefilteringResult:
        if both_strains:
            return self.prefiltering.calculate_top_matches_with_rev_comp(query)
        return self.prefiltering.calculate_top_matches(query)

    def _align_sequences(self, query: str, prefiltering_result: SpeciesPrefilteringResult) -> list[AlignmentEntryNT]:
        aligner = StripedSmithWaterman(query, **ALIGNER_PARAMS)
        rev_comp_aligner = StripedSmithWaterman(prefiltering_result.rev_comp_query, **ALIGNER_PARAMS)

        result = []

        for gene_match in prefiltering_result.top_matches:
            current_species = gene_match.species
            match gene_match.rev_comp:
                case False:
                    alignment = align(
                        aligner=aligner,
                        target_id=gene_match.gene_id,
                        target=self.gene_lookup[gene_match.species + "|" + gene_match.gene_id].sequence,
                        query=query,
                        db_length=self.db_length,
                        rev_comp=gene_match.rev_comp,
                    )
                case True:
                    alignment = align(
                        aligner=rev_comp_aligner,
                        target_id=gene_match.gene_id,
                        target=self.gene_lookup[gene_match.species + "|" + gene_match.gene_id].sequence,
                        query=prefiltering_result.rev_comp_query,
                        db_length=self.db_length,
                        rev_comp=gene_match.rev_comp,
                    )

            alignment_length = alignment.t_end - alignment.t_start

            if self.e_value_threshold is not None:
                assert alignment.e_value is not None
                if alignment.e_value > self.e_value_threshold:
                    continue
            if self.alignment_length_threshold is not None:
                if alignment_length < self.alignment_length_threshold:
                    continue

            target_gene: Gene = self.gene_lookup[gene_match.species + "|" + gene_match.gene_id]

            result.append(
                AlignmentEntryNT(
                    target_id=alignment.target_id,
                    alignment_score=alignment.alignment_score,
                    seq_identity=alignment.seq_identity,
                    e_value=alignment.e_value,
                    q_start=alignment.q_start,
                    q_end=alignment.q_end,
                    q_len=alignment.q_len,
                    t_start=alignment.t_start,
                    t_end=alignment.t_end,
                    t_len=alignment.t_len,
                    cigar=alignment.cigar,
                    rev_comp=alignment.rev_comp,
                    species=current_species,
                    locus=gene_match.locus,
                    q_seq=alignment.query,
                    t_seq=target_gene.sequence,
                    reading_frame=target_gene.reading_frame,
                )
            )

        return result

    def align(self, query: str, both_strains: bool) -> Optional[AlignmentEntryNT]:
        prefiltering_result = self._prefilter(query, both_strains)
        alignments = self._align_sequences(query, prefiltering_result)

        alignments.sort()

        return alignments[0] if len(alignments) > 0 else None


AA_ALIGNER_PARAMS = {
    "gap_open_penalty": 11,
    "gap_extend_penalty": 1,
    "protein": True,
    "substitution_matrix": blosum.BLOSUM(62),
}


class GeneAlignerAA:
    def __init__(
        self,
        genes: list[Gene],
        kmer_size: int,
        distance_threshold: int,
        top_n: int,
        modulo_n: int,
        e_value_threshold: Optional[float] = None,
        alignment_length_threshold: Optional[int] = None,
        max_cdr3_length: Optional[int] = None,
    ) -> None:
        self.prefiltering = MultiSpeciesPrefiltering(genes, kmer_size, distance_threshold, top_n, modulo_n)

        self.gene_lookup: dict[str, Gene] = dict(map(lambda gene: (gene.species + "|" + gene.name, gene), genes))
        self.db_length = sum(len(gene.sequence) for gene in genes)

        self.e_value_threshold = e_value_threshold
        self.alignment_length_threshold = alignment_length_threshold
        self.max_cdr3_length = max_cdr3_length

    def _prefilter(self, query: str) -> SpeciesPrefilteringResult:
        return self.prefiltering.calculate_top_matches(query)

    def _align_sequences(self, query: str, prefiltering_result: SpeciesPrefilteringResult) -> list[AlignmentEntryAA]:
        aligner = StripedSmithWaterman(query, **AA_ALIGNER_PARAMS)

        result = []

        for gene_match in prefiltering_result.top_matches:
            current_species = gene_match.species
            alignment = align_aa(
                aligner=aligner,
                target_id=gene_match.gene_id,
                target=self.gene_lookup[gene_match.species + "|" + gene_match.gene_id].sequence,
                query=query,
                db_length=self.db_length,
            )

            alignment_length = alignment.t_end - alignment.t_start

            if self.e_value_threshold is not None:
                assert alignment.e_value is not None
                if alignment.e_value > self.e_value_threshold:
                    continue
            if self.alignment_length_threshold is not None:
                if alignment_length < self.alignment_length_threshold:
                    continue
            if self.max_cdr3_length is not None:
                if alignment.q_start > self.max_cdr3_length:
                    continue

            target_gene: Gene = self.gene_lookup[gene_match.species + "|" + gene_match.gene_id]

            result.append(
                AlignmentEntryAA(
                    target_id=alignment.target_id,
                    alignment_score=alignment.alignment_score,
                    seq_identity=alignment.seq_identity,
                    e_value=alignment.e_value,
                    q_start=alignment.q_start,
                    q_end=alignment.q_end,
                    t_start=alignment.t_start,
                    t_end=alignment.t_end,
                    cigar=alignment.cigar,
                    species=current_species,
                    locus=gene_match.locus,
                    q_seq=query,
                    t_seq=target_gene.sequence,
                )
            )

        return result

    def align(self, query: str) -> Optional[AlignmentEntryAA]:
        prefiltering_result = self._prefilter(query)
        alignments = self._align_sequences(query, prefiltering_result)

        alignments.sort()

        return alignments[0] if len(alignments) > 0 else None


def get_aligner_params(germline_gene: GermlineGene, locus: Optional[Locus]) -> dict:
    # ideally this should be read from a file or env

    match germline_gene:
        case GermlineGene.V:
            return {
                "kmer_size": 9,
                "distance_threshold": 13,
                "top_n": int(os.environ.get("TOP_N", 12)),
                "modulo_n": 2,
                "e_value_threshold": 0.001,
                "alignment_length_threshold": 100,
            }
        case GermlineGene.D:
            return {
                "kmer_size": 5,
                "distance_threshold": 3,
                "top_n": int(os.environ.get("TOP_N", 5)),
                "modulo_n": 1,
            }
        case GermlineGene.J:
            assert locus is not None
            match locus:
                case Locus.IGH:
                    return {
                        "kmer_size": 5,
                        "distance_threshold": 3,
                        "top_n": int(os.environ.get("TOP_N", 5)),
                        "modulo_n": 2,
                    }
                case Locus.IGK:
                    return {
                        "kmer_size": 5,
                        "distance_threshold": 5,
                        "top_n": int(os.environ.get("TOP_N", 5)),
                        "modulo_n": 1,
                    }
                case Locus.IGL:
                    return {
                        "kmer_size": 5,
                        "distance_threshold": 7,
                        "top_n": int(os.environ.get("TOP_N", 5)),
                        "modulo_n": 2,
                    }
        case GermlineGene.C:
            assert locus is not None
            match locus:
                case Locus.IGH:
                    return {
                        "kmer_size": 5,
                        "distance_threshold": 5,
                        "top_n": int(os.environ.get("TOP_N", 5)),
                        "modulo_n": 2,
                    }
                case Locus.IGK:
                    return {
                        "kmer_size": 3,
                        "distance_threshold": 5,
                        "top_n": int(os.environ.get("TOP_N", 3)),
                        "modulo_n": 1,
                    }
                case Locus.IGL:
                    return {
                        "kmer_size": 3,
                        "distance_threshold": 5,
                        "top_n": int(os.environ.get("TOP_N", 3)),
                        "modulo_n": 1,
                    }


def get_aa_aligner_params(germline_gene: GermlineGene) -> dict:
    # ideally this should be read from a file or env

    match germline_gene:
        case GermlineGene.V:
            return {
                "kmer_size": 3,
                "distance_threshold": 3,
                "top_n": int(os.environ.get("TOP_N", 20)),
                "modulo_n": 1,
                "e_value_threshold": 1e-55,
                "alignment_length_threshold": 80,
            }
        case GermlineGene.J:
            return {
                "kmer_size": 3,
                "distance_threshold": 3,
                "top_n": int(os.environ.get("TOP_N", 5)),
                "modulo_n": 1,
                "alignment_length_threshold": 7,
                "max_cdr3_length": 60,
            }
        case _:
            raise ValueError("AA aligner not supported for gene type")


def get_gene_parsing_function(gene_type: GermlineGene) -> Callable[[str, SeqRecord], Gene]:
    match gene_type:
        case GermlineGene.V:
            return lambda description, record: Gene(
                name=description[0],
                locus=Locus[description[1]],
                reading_frame=int(description[2]),
                species=Organism[description[3]],
                sequence=str(record.seq),
            )
        case GermlineGene.D:
            return lambda description, record: Gene(
                name=description[0],
                locus=Locus.IGH,
                reading_frame=None,
                species=Organism[description[1]],
                sequence=str(record.seq),
            )
        case GermlineGene.J:
            return lambda description, record: Gene(
                name=description[0],
                locus=Locus[description[1]],
                reading_frame=int(description[2]),
                species=Organism[description[3]],
                sequence=str(record.seq),
            )
        case GermlineGene.C:
            return lambda description, record: Gene(
                name=description[0],
                locus=Locus[description[1]],
                reading_frame=None,
                species=Organism.HOMO_SAPIENS,
                sequence=str(record.seq),
            )


def parse_gene_aa(description: Sequence[str], record: SeqRecord, species: Organism):
    return GeneAA(name=description[0], locus=Locus[description[1]], species=species, sequence=str(record.seq))


def read_genes(path: Path, parsing_fn):
    result = []
    for record in SeqIO.parse(path, "fasta"):
        description = record.description.split("\t")
        gene = parsing_fn(description, record)

        result.append(gene)

    return result


def create_v_gene_aligner(allowed_species: list[Organism], db_dir: Path = GENE_DB_DIR) -> GeneAligner:
    genes = []

    if not allowed_species:
        allowed_species = [Organism.HOMO_SAPIENS, Organism.MUS_MUSCULUS]

    gene_parsing_function = get_gene_parsing_function(GermlineGene.V)

    for species in allowed_species:
        input_path = Path(db_dir) / "gene_db" / "v_genes" / f"{species.value}.fasta"
        species_genes = read_genes(input_path, gene_parsing_function)

        genes.extend(species_genes)

    aligner_params = get_aligner_params(GermlineGene.V, None)
    genes_aligner = GeneAligner(genes, **aligner_params)

    return genes_aligner


def create_aligner(
    allowed_species: Organism, germline_gene: GermlineGene, locus: Locus, db_dir: Path = GENE_DB_DIR
) -> GeneAligner:
    assert germline_gene != GermlineGene.V

    gene_parsing_function = get_gene_parsing_function(germline_gene)

    gene_path = (
        Path(db_dir)
        / "gene_db"
        / f"{germline_gene.value.lower()}_genes"
        / allowed_species.value
        / f"{locus.value}.fasta"
    )
    genes = read_genes(gene_path, gene_parsing_function)

    aligner_params = get_aligner_params(germline_gene, locus)
    genes_aligner = GeneAligner(genes, **aligner_params)

    return genes_aligner


def create_aa_v_gene_aligner(allowed_species: list[Organism], aa_genes_dir: Path) -> GeneAlignerAA:
    genes = []

    if not allowed_species:
        allowed_species = [Organism.HOMO_SAPIENS, Organism.MUS_MUSCULUS]

    for species in allowed_species:
        input_path = aa_genes_dir / "v_genes" / f"{species.value}.fasta"

        assert input_path.exists(), f"Input germline file: {str(input_path)} does not exists"

        species_genes = read_genes(input_path, partial(parse_gene_aa, species=species))

        genes.extend(species_genes)

    aligner_params = get_aa_aligner_params(GermlineGene.V)
    genes_aligner = GeneAlignerAA(genes, **aligner_params)

    return genes_aligner


def create_aa_j_gene_aligner(organism: Organism, locus: Locus, aa_genes_dir: Path) -> GeneAlignerAA:
    input_path = aa_genes_dir / "j_genes" / organism.value / f"{locus.value}.fasta"
    assert input_path.exists(), f"Input germline file: {str(input_path)} does not exists"

    genes = read_genes(input_path, partial(parse_gene_aa, species=organism))

    aligner_params = get_aa_aligner_params(GermlineGene.J)
    genes_aligner = GeneAlignerAA(genes, **aligner_params)

    return genes_aligner
