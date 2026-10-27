import math
import os
from functools import cached_property, partial
from pathlib import Path
from typing import Callable, Optional, Sequence

from ssw_aligner import StripedSmithWaterman

from riot_na.alignment.alignment_metrics import GumbellParams, compute_bit_score
from riot_na.alignment.skbio_alignment import (
    align_aa,  # type: ignore  # pylint: disable=import-error
)
from riot_na.alignment.skbio_alignment import align
from riot_na.common.fasta import FastaRecord, parse_fasta
from riot_na.common.gene_match_utils import create_gene_lookup, gene_lookup_key
from riot_na.common.multi_species_prefiltering import MultiSpeciesPrefiltering
from riot_na.config import GENE_DB_DIR
from riot_na.data.constants import BLOSUM_62
from riot_na.data.model import (
    AlignmentEntryAA,
    AlignmentEntryNT,
    Gene,
    GeneAA,
    GermlineGene,
    Locus,
    Organism,
    SpeciesGeneMatch,
    SpeciesPrefilteringResult,
)

# germlines selected by name are looked up in human germlines first, if more species have germlines of the same name
GENE_SELECTION_SPECIES_ORDER = (
    Organism.HOMO_SAPIENS,
    *(species for species in Organism if species != Organism.HOMO_SAPIENS),
)

ALIGNER_PARAMS = {"match_score": 1, "mismatch_score": -1, "gap_open_penalty": 4, "gap_extend_penalty": 1}


class GeneAligner:
    def __init__(
        self,
        genes: Sequence[Gene],
        kmer_size: int,
        distance_threshold: int,
        top_n: int,
        modulo_n: int,
        e_value_threshold: Optional[float] = None,
        alignment_length_threshold: Optional[int] = None,
    ) -> None:

        self.prefiltering = MultiSpeciesPrefiltering(genes, kmer_size, distance_threshold, top_n, modulo_n)

        self.gene_lookup = create_gene_lookup(genes)
        self.db_length = sum(len(gene.sequence) for gene in genes)

        self.e_value_threshold = e_value_threshold
        self.alignment_length_threshold = alignment_length_threshold

    def _prefilter(self, query: str, both_strains: bool) -> SpeciesPrefilteringResult:
        """Prefilter genes and return species-aware gene matches."""

        if both_strains:
            return self.prefiltering.calculate_top_matches_with_rev_comp(query)

        return self.prefiltering.calculate_top_matches(query)

    def _align_sequences(self, query: str, prefiltering_result: SpeciesPrefilteringResult) -> list[AlignmentEntryNT]:
        aligner = StripedSmithWaterman(query, **ALIGNER_PARAMS)
        rev_comp_aligner = (
            StripedSmithWaterman(prefiltering_result.rev_comp_query, **ALIGNER_PARAMS)
            if prefiltering_result.rev_comp_query
            else None
        )

        result = []

        for gene_match in prefiltering_result.top_matches:
            target_gene = self.gene_lookup[gene_match.species_gene_id]
            current_query = (
                prefiltering_result.rev_comp_query
                if gene_match.rev_comp and prefiltering_result.rev_comp_query
                else query
            )
            current_aligner = rev_comp_aligner if gene_match.rev_comp and rev_comp_aligner else aligner

            alignment = align(
                aligner=current_aligner,
                target_id=gene_match.gene_id,
                target=target_gene.sequence,
                query=current_query,
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
                    species=gene_match.species,
                    locus=gene_match.locus,
                    q_seq=current_query,
                    t_seq=target_gene.sequence,
                    reading_frame=target_gene.reading_frame if isinstance(target_gene, Gene) else None,
                )
            )

        return result

    def align(self, query: str, both_strains: bool) -> Optional[AlignmentEntryNT]:
        """Best alignment (lowest e-value) of the query to the genes."""
        prefiltering_result = self._prefilter(query, both_strains)
        alignments = self._align_sequences(query, prefiltering_result)
        return min(alignments, default=None)

    def max_alignment_score(self, query: str) -> float:
        """
        Upper bound of alignment scores (bit scores) of the query to the genes: aligned nucleotides score at most
        match_score, mismatches and gaps lower the score.
        """
        return compute_bit_score(ALIGNER_PARAMS["match_score"] * min(len(query), self._max_gene_length))

    @cached_property
    def _max_gene_length(self) -> int:
        return max((len(gene.sequence) for gene in self.gene_lookup.values()), default=0)


AA_ALIGNER_PARAMS = {
    "gap_open_penalty": 11,
    "gap_extend_penalty": 1,
    "protein": True,
    "substitution_matrix": BLOSUM_62,
}

# best substitution score of every residue, residues missing in the matrix get the best score of the matrix
_MAX_SUBSTITUTION_SCORES_AA = {
    residue: max(0, math.ceil(max(scores.values()))) for residue, scores in BLOSUM_62.items()
}
_MAX_SUBSTITUTION_SCORE_AA = max(_MAX_SUBSTITUTION_SCORES_AA.values())


def max_raw_alignment_score_aa(sequence: str) -> int:
    """
    Upper bound of raw scores of local alignments of the sequence (to any sequence): aligned residues score at most
    their best substitution score, gaps lower the score.
    """
    return sum(_MAX_SUBSTITUTION_SCORES_AA.get(residue, _MAX_SUBSTITUTION_SCORE_AA) for residue in sequence)


class GeneAlignerAA:  # pylint: disable=too-many-instance-attributes
    def __init__(  # pylint: disable=too-many-arguments
        self,
        genes: Sequence[GeneAA],
        kmer_size: int,
        distance_threshold: int,
        top_n: int,
        modulo_n: int,
        e_value_threshold: Optional[float] = None,
        alignment_length_threshold: Optional[int] = None,
        max_cdr3_length: Optional[int] = None,
        min_seq_identity: Optional[float] = None,
        min_internal_alignment_length: Optional[int] = None,
        max_fragment_germline_start: int = 0,
    ) -> None:
        # Initialize the base prefiltering directly with species-aware keys
        self.prefiltering = MultiSpeciesPrefiltering(genes, kmer_size, distance_threshold, top_n, modulo_n)

        self.gene_lookup = create_gene_lookup(genes)
        self.db_length = sum(len(gene.sequence) for gene in genes)

        self.e_value_threshold = e_value_threshold
        self.alignment_length_threshold = alignment_length_threshold
        self.max_cdr3_length = max_cdr3_length
        self.min_seq_identity = min_seq_identity
        self.min_internal_alignment_length = min_internal_alignment_length
        self.max_fragment_germline_start = max_fragment_germline_start

    def find_gene(self, gene_id: str) -> Optional[Gene | GeneAA]:
        """
        Look up a germline by its exact name, e.g. IGHV3-7*01, in the species the aligner was created for.
        Human germline is returned if more species have a germline of this name.
        """
        for species in GENE_SELECTION_SPECIES_ORDER:
            for locus in Locus:
                gene = self.gene_lookup.get(gene_lookup_key(species, locus, gene_id))
                if gene is not None:
                    return gene
        return None

    def get_gene(self, gene_id: str) -> Gene | GeneAA:
        gene = self.find_gene(gene_id)
        assert gene is not None, f"Gene {gene_id} not found in germline database"
        return gene

    def _prefilter(self, query: str) -> SpeciesPrefilteringResult:
        """Prefilter genes and return species-aware gene matches."""
        return self.prefiltering.calculate_top_matches(query)

    def _align_sequences(
        self, query: str, prefiltering_result: SpeciesPrefilteringResult
    ) -> Sequence[AlignmentEntryAA]:
        aligner = StripedSmithWaterman(query, **AA_ALIGNER_PARAMS)

        result = []

        for gene_match in prefiltering_result.top_matches:
            target_gene = self.gene_lookup[gene_match.species_gene_id]

            alignment = align_aa(
                aligner=aligner,
                target_id=gene_match.gene_id,
                target=target_gene.sequence,
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
            if self.min_seq_identity is not None:
                assert alignment.seq_identity is not None
                if alignment.seq_identity < self.min_seq_identity:
                    continue
            if self.min_internal_alignment_length is not None:
                # alignments shorter than that have to start at the beginning of the germline, or at most
                # max_fragment_germline_start residues after it if they start at the beginning of the query too
                if alignment_length < self.min_internal_alignment_length and (
                    alignment.t_start > self.max_fragment_germline_start
                    or (alignment.t_start > 0 and alignment.q_start > self.max_fragment_germline_start)
                ):
                    continue

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
                    species=gene_match.species,
                    locus=gene_match.locus,
                    q_seq=query,
                    t_seq=target_gene.sequence,
                )
            )

        return result

    def align(self, query: str) -> Optional[AlignmentEntryAA]:
        """Best alignment (lowest e-value) of the query to the genes."""
        prefiltering_result = self._prefilter(query)
        alignments = self._align_sequences(query, prefiltering_result)
        return min(alignments, default=None)

    def max_alignment_score(self, query: str) -> float:
        """Upper bound of alignment scores (bit scores) of the query to the genes (see max_raw_alignment_score_aa)."""
        return compute_bit_score(
            min(self._max_gene_raw_score, max_raw_alignment_score_aa(query)), gumbell_params=GumbellParams.AA
        )

    @cached_property
    def _max_gene_raw_score(self) -> int:
        return max((max_raw_alignment_score_aa(gene.sequence) for gene in self.gene_lookup.values()), default=0)

    def align_to_genes(self, query: str, gene_ids: Sequence[str]) -> Sequence[AlignmentEntryAA]:
        """
        Align query to explicitly selected germlines, skipping prefiltering.

        Alignment thresholds are applied just like in the regular, prefiltered path,
        so a poorly matching germline may yield no alignment at all.
        """
        gene_matches = []

        for gene_id in gene_ids:
            gene = self.get_gene(gene_id)
            gene_matches.append(
                SpeciesGeneMatch(
                    species_gene_id=gene_lookup_key(gene.species, gene.locus, gene.name),
                    gene_id=gene.name,
                    rev_comp=False,
                    coverage=0,
                    species=gene.species,
                    locus=gene.locus,
                )
            )

        prefiltering_result = SpeciesPrefilteringResult(query=query, rev_comp_query="", top_matches=gene_matches)

        return self._align_sequences(query, prefiltering_result)


def get_aligner_params(germline_gene: GermlineGene, locus: Optional[Locus]) -> dict:
    # ideally this should be read from a file or env
    match germline_gene:
        case GermlineGene.V:
            return {
                "kmer_size": int(os.environ.get("KMER_SIZE_V_NT", 5)),
                "distance_threshold": int(os.environ.get("DISTANCE_THRESHOLD_V_NT", 13)),
                "top_n": int(os.environ.get("TOP_N_V_NT", 12)),
                "modulo_n": int(os.environ.get("MODULO_N_V_NT", 2)),
                "e_value_threshold": float(os.environ.get("E_VALUE_THRESHOLD_V_NT", 0.05)),
                "alignment_length_threshold": int(os.environ.get("ALIGNMENT_LENGTH_THRESHOLD_V_NT", 100)),
            }
        case GermlineGene.D:
            return {
                "kmer_size": int(os.environ.get("KMER_SIZE_D_NT", 5)),
                "distance_threshold": int(os.environ.get("DISTANCE_THRESHOLD_D_NT", 3)),
                "top_n": int(os.environ.get("TOP_N_D_NT", 5)),
                "modulo_n": int(os.environ.get("MODULO_N_D_NT", 1)),
            }
        case GermlineGene.J:
            assert locus is not None
            match locus:
                case Locus.IGH:
                    return {
                        "kmer_size": int(os.environ.get("KMER_SIZE_JH_NT", 5)),
                        "distance_threshold": int(os.environ.get("DISTANCE_THRESHOLD_JH_NT", 3)),
                        "top_n": int(os.environ.get("TOP_N_JH_NT", 5)),
                        "modulo_n": int(os.environ.get("MODULO_N_JH_NT", 2)),
                    }
                case Locus.IGK:
                    return {
                        "kmer_size": int(os.environ.get("KMER_SIZE_JK_NT", 5)),
                        "distance_threshold": int(os.environ.get("DISTANCE_THRESHOLD_JK_NT", 5)),
                        "top_n": int(os.environ.get("TOP_N_JK_NT", 5)),
                        "modulo_n": int(os.environ.get("MODULO_N_JK_NT", 1)),
                    }
                case Locus.IGL:
                    return {
                        "kmer_size": int(os.environ.get("KMER_SIZE_JL_NT", 5)),
                        "distance_threshold": int(os.environ.get("DISTANCE_THRESHOLD_JL_NT", 7)),
                        "top_n": int(os.environ.get("TOP_N_JL_NT", 5)),
                        "modulo_n": int(os.environ.get("MODULO_N_JL_NT", 2)),
                    }
        case GermlineGene.C:
            assert locus is not None
            match locus:
                case Locus.IGH:
                    return {
                        "kmer_size": int(os.environ.get("KMER_SIZE_CH_NT", 5)),
                        "distance_threshold": int(os.environ.get("DISTANCE_THRESHOLD_CH_NT", 5)),
                        "top_n": int(os.environ.get("TOP_N_CH_NT", 5)),
                        "modulo_n": int(os.environ.get("MODULO_N_CH_NT", 2)),
                        "e_value_threshold": float(os.environ.get("E_VALUE_THRESHOLD_CH_NT", 0.05)),
                    }
                case Locus.IGK:
                    return {
                        "kmer_size": int(os.environ.get("KMER_SIZE_CK_NT", 3)),
                        "distance_threshold": int(os.environ.get("DISTANCE_THRESHOLD_CK_NT", 5)),
                        "top_n": int(os.environ.get("TOP_N_CK_NT", 3)),
                        "modulo_n": int(os.environ.get("MODULO_N_CK_NT", 1)),
                        "e_value_threshold": float(os.environ.get("E_VALUE_THRESHOLD_CK_NT", 0.05)),
                    }
                case Locus.IGL:
                    return {
                        "kmer_size": int(os.environ.get("KMER_SIZE_CL_NT", 3)),
                        "distance_threshold": int(os.environ.get("DISTANCE_THRESHOLD_CL_NT", 5)),
                        "top_n": int(os.environ.get("TOP_N_CL_NT", 3)),
                        "modulo_n": int(os.environ.get("MODULO_N_CL_NT", 1)),
                        "e_value_threshold": float(os.environ.get("E_VALUE_THRESHOLD_CL_NT", 0.05)),
                    }


def get_aa_aligner_params(germline_gene: GermlineGene) -> dict:
    # ideally this should be read from a file or env

    match germline_gene:
        case GermlineGene.V:
            return {
                "kmer_size": int(os.environ.get("KMER_SIZE_V_AA", 3)),
                "distance_threshold": int(os.environ.get("DISTANCE_THRESHOLD_V_AA", 3)),
                "top_n": int(os.environ.get("TOP_N_V_AA", 20)),
                "modulo_n": int(os.environ.get("MODULO_N_V_AA", 1)),
                "e_value_threshold": float(os.environ.get("E_VALUE_THRESHOLD_V_AA", 1e-45)),
                "alignment_length_threshold": int(os.environ.get("ALIGNMENT_LENGTH_THRESHOLD_V_AA", 50)),
            }
        case GermlineGene.J:
            return {
                "kmer_size": int(os.environ.get("KMER_SIZE_J_AA", 3)),
                "distance_threshold": int(os.environ.get("DISTANCE_THRESHOLD_J_AA", 3)),
                "top_n": int(os.environ.get("TOP_N_J_AA", 5)),
                "modulo_n": int(os.environ.get("MODULO_N_J_AA", 1)),
                "alignment_length_threshold": int(os.environ.get("ALIGNMENT_LENGTH_THRESHOLD_J_AA", 5)),
                "e_value_threshold": float(os.environ.get("E_VALUE_THRESHOLD_J_AA", 0.001)),
                "max_cdr3_length": int(os.environ.get("MAX_CDR3_LENGTH_J_AA", 60)),
            }
        case GermlineGene.C:
            return {
                "kmer_size": int(os.environ.get("KMER_SIZE_C_AA", 3)),
                "distance_threshold": int(os.environ.get("DISTANCE_THRESHOLD_C_AA", 3)),
                "top_n": int(os.environ.get("TOP_N_C_AA", 5)),
                "modulo_n": int(os.environ.get("MODULO_N_C_AA", 2)),
                "alignment_length_threshold": int(os.environ.get("ALIGNMENT_LENGTH_THRESHOLD_C_AA", 4)),
                "e_value_threshold": float(os.environ.get("E_VALUE_THRESHOLD_C_AA", 0.05)),
                # filter out matches of non antibody regions following V region, e.g. proline rich CAR hinges
                # to IgG and IgA hinges. Short C regions (e.g. CH1 fragments in Fabs) start at the C germline start.
                "min_seq_identity": float(os.environ.get("MIN_SEQ_IDENTITY_C_AA", 0.4)),
                "min_internal_alignment_length": int(os.environ.get("MIN_INTERNAL_ALIGNMENT_LENGTH_C_AA", 50)),
                # the first C residue is encoded by the J-C junction codon (its first nucleotide comes from J), so it is
                # missing or different when J ends differently, e.g. absorbed by J alignment or mutated. C fragments
                # without it start right after J (or at the start of the query of the constant region aligner)
                "max_fragment_germline_start": int(os.environ.get("MAX_FRAGMENT_GERMLINE_START_C_AA", 1)),
            }
        case _:
            raise ValueError("AA aligner not supported for gene type")


def get_gene_parsing_function(gene_type: GermlineGene) -> Callable[[Sequence[str], FastaRecord], Gene]:
    match gene_type:
        case GermlineGene.V:
            return lambda description, record: Gene(
                name=description[0],
                locus=Locus[description[1]],
                reading_frame=int(description[2]),
                species=Organism[description[3]],
                sequence=record.seq,
            )
        case GermlineGene.D:
            return lambda description, record: Gene(
                name=description[0],
                locus=Locus.IGH,
                reading_frame=None,
                species=Organism[description[1]],
                sequence=record.seq,
            )
        case GermlineGene.J:
            return lambda description, record: Gene(
                name=description[0],
                locus=Locus[description[1]],
                reading_frame=int(description[2]),
                species=Organism[description[3]],
                sequence=record.seq,
            )
        case GermlineGene.C:
            return lambda description, record: Gene(
                name=description[0],
                locus=Locus[description[1]],
                reading_frame=None,
                species=Organism[description[2]],
                sequence=record.seq,
            )


def parse_gene_aa(description: Sequence[str], record: FastaRecord, species: Organism):
    return GeneAA(name=description[0], locus=Locus[description[1]], species=species, sequence=record.seq)


def read_genes(path: Path, parsing_fn):
    result = []
    for record in parse_fasta(path):
        description = record.description.split("\t")
        gene = parsing_fn(description, record)

        result.append(gene)

    return result


def create_v_gene_aligner(allowed_species: Sequence[Organism], db_dir: Path = GENE_DB_DIR) -> GeneAligner:
    genes = []

    if not allowed_species:
        allowed_species = (Organism.HOMO_SAPIENS, Organism.MUS_MUSCULUS, Organism.VICUGNA_PACOS)

    gene_parsing_function = get_gene_parsing_function(GermlineGene.V)

    for species in allowed_species:
        input_path = Path(db_dir) / "gene_db" / "v_genes" / f"{species.value}.fasta"
        species_genes = read_genes(input_path, gene_parsing_function)

        genes.extend(species_genes)

    aligner_params = get_aligner_params(GermlineGene.V, None)
    genes_aligner = GeneAligner(genes, **aligner_params)

    return genes_aligner


def create_aligner(
    organism: Organism, germline_gene: GermlineGene, locus: Locus, db_dir: Path = GENE_DB_DIR
) -> GeneAligner:
    assert germline_gene != GermlineGene.V

    gene_parsing_function = get_gene_parsing_function(germline_gene)

    gene_path = (
        Path(db_dir) / "gene_db" / f"{germline_gene.value.lower()}_genes" / organism.value / f"{locus.value}.fasta"
    )
    genes = read_genes(gene_path, gene_parsing_function)

    aligner_params = get_aligner_params(germline_gene, locus)
    genes_aligner = GeneAligner(genes, **aligner_params)

    return genes_aligner


def create_aa_v_gene_aligner(allowed_species: Sequence[Organism], aa_genes_dir: Path) -> GeneAlignerAA:
    genes = []

    if not allowed_species:
        allowed_species = (Organism.HOMO_SAPIENS, Organism.MUS_MUSCULUS, Organism.VICUGNA_PACOS)

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


def create_aa_c_gene_aligner(organism: Organism, locus: Locus, aa_genes_dir: Path) -> GeneAlignerAA:
    input_path = aa_genes_dir / "c_genes" / organism.value / f"{locus.value}.fasta"
    assert input_path.exists(), f"Input germline file: {str(input_path)} does not exists"

    genes = read_genes(input_path, partial(parse_gene_aa, species=organism))

    aligner_params = get_aa_aligner_params(GermlineGene.C)
    genes_aligner = GeneAlignerAA(genes, **aligner_params)

    return genes_aligner
