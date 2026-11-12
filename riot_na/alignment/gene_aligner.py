import math
from functools import cached_property, partial
from pathlib import Path
from typing import Callable, Optional, Sequence

from ssw_aligner import StripedSmithWaterman

from riot_na.alignment.alignment_metrics import (
    compute_bit_score,
    log10_evalue_threshold,
)
from riot_na.alignment.gumbel import get_gumbel_params
from riot_na.alignment.params import (
    DEFAULT_ALIGNER_PARAMS,
    AASSWParams,
    AlignerParams,
    NTSSWParams,
)
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


class GeneAligner:  # pylint: disable=too-many-instance-attributes
    def __init__(
        self,
        genes: Sequence[Gene],
        kmer_size: int,
        distance_threshold: int,
        top_n: int,
        modulo_n: int,
        e_value_threshold: Optional[float] = None,
        alignment_length_threshold: Optional[int] = None,
        min_seq_identity: Optional[float] = None,
        ssw_params: NTSSWParams = DEFAULT_ALIGNER_PARAMS.v_nt.aligner,
    ) -> None:

        self.ssw_params = ssw_params
        self.gumbel_params = get_gumbel_params(ssw_params)
        self.prefiltering = MultiSpeciesPrefiltering(genes, kmer_size, distance_threshold, top_n, modulo_n)

        self.gene_lookup = create_gene_lookup(genes)
        self.db_length = sum(len(gene.sequence) for gene in genes)

        self.e_value_threshold = e_value_threshold
        self.log10_e_value_threshold = log10_evalue_threshold(e_value_threshold)
        self.alignment_length_threshold = alignment_length_threshold
        self.min_seq_identity = min_seq_identity

    def _prefilter(self, query: str, both_strains: bool) -> SpeciesPrefilteringResult:
        """Prefilter genes and return species-aware gene matches."""

        if both_strains:
            return self.prefiltering.calculate_top_matches_with_rev_comp(query)

        return self.prefiltering.calculate_top_matches(query)

    def _align_sequences(self, query: str, prefiltering_result: SpeciesPrefilteringResult) -> list[AlignmentEntryNT]:
        aligner = StripedSmithWaterman(query, **self.ssw_params.kwargs())
        rev_comp_aligner = (
            StripedSmithWaterman(prefiltering_result.rev_comp_query, **self.ssw_params.kwargs())
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
                gumbel_params=self.gumbel_params,
                rev_comp=gene_match.rev_comp,
            )

            alignment_length = alignment.t_end - alignment.t_start

            if self.log10_e_value_threshold is not None:
                assert alignment.log10_e_value is not None
                if alignment.log10_e_value > self.log10_e_value_threshold:
                    continue
            if self.alignment_length_threshold is not None:
                if alignment_length < self.alignment_length_threshold:
                    continue
            if self.min_seq_identity is not None:
                assert alignment.seq_identity is not None
                if alignment.seq_identity < self.min_seq_identity:
                    continue

            result.append(
                AlignmentEntryNT(
                    target_id=alignment.target_id,
                    alignment_score=alignment.alignment_score,
                    seq_identity=alignment.seq_identity,
                    log10_e_value=alignment.log10_e_value,
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
        return compute_bit_score(
            self.ssw_params.match_score * min(len(query), self._max_gene_length), gumbel_params=self.gumbel_params
        )

    @cached_property
    def _max_gene_length(self) -> int:
        return max((len(gene.sequence) for gene in self.gene_lookup.values()), default=0)


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
        ssw_params: AASSWParams = DEFAULT_ALIGNER_PARAMS.v_aa.aligner,
    ) -> None:
        self.ssw_params = ssw_params
        self.gumbel_params = get_gumbel_params(ssw_params)
        # Initialize the base prefiltering directly with species-aware keys
        self.prefiltering = MultiSpeciesPrefiltering(genes, kmer_size, distance_threshold, top_n, modulo_n)

        self.gene_lookup = create_gene_lookup(genes)
        self.db_length = sum(len(gene.sequence) for gene in genes)

        self.e_value_threshold = e_value_threshold
        self.log10_e_value_threshold = log10_evalue_threshold(e_value_threshold)
        self.alignment_length_threshold = alignment_length_threshold
        self.max_cdr3_length = max_cdr3_length
        self.min_seq_identity = min_seq_identity
        self.min_internal_alignment_length = min_internal_alignment_length

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
        aligner = StripedSmithWaterman(query, **self.ssw_params.kwargs())

        result = []

        for gene_match in prefiltering_result.top_matches:
            target_gene = self.gene_lookup[gene_match.species_gene_id]

            alignment = align_aa(
                aligner=aligner,
                target_id=gene_match.gene_id,
                target=target_gene.sequence,
                query=query,
                db_length=self.db_length,
                gumbel_params=self.gumbel_params,
            )

            alignment_length = alignment.t_end - alignment.t_start

            if self.log10_e_value_threshold is not None:
                assert alignment.log10_e_value is not None
                if alignment.log10_e_value > self.log10_e_value_threshold:
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
                if alignment_length < self.min_internal_alignment_length:
                    continue

            result.append(
                AlignmentEntryAA(
                    target_id=alignment.target_id,
                    alignment_score=alignment.alignment_score,
                    seq_identity=alignment.seq_identity,
                    log10_e_value=alignment.log10_e_value,
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
            min(self._max_gene_raw_score, max_raw_alignment_score_aa(query)), gumbel_params=self.gumbel_params
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


def create_v_gene_aligner(
    allowed_species: Sequence[Organism],
    db_dir: Path = GENE_DB_DIR,
    params: AlignerParams = DEFAULT_ALIGNER_PARAMS,
) -> GeneAligner:
    genes = []

    if not allowed_species:
        allowed_species = (Organism.HOMO_SAPIENS, Organism.MUS_MUSCULUS, Organism.VICUGNA_PACOS)

    gene_parsing_function = get_gene_parsing_function(GermlineGene.V)

    for species in allowed_species:
        input_path = Path(db_dir) / "gene_db" / "v_genes" / f"{species.value}.fasta"
        species_genes = read_genes(input_path, gene_parsing_function)

        genes.extend(species_genes)

    genes_aligner = GeneAligner(genes, **params.v_nt.kwargs(), ssw_params=params.v_nt.aligner)

    return genes_aligner


def create_aligner(
    organism: Organism,
    germline_gene: GermlineGene,
    locus: Locus,
    db_dir: Path = GENE_DB_DIR,
    params: AlignerParams = DEFAULT_ALIGNER_PARAMS,
) -> GeneAligner:
    assert germline_gene != GermlineGene.V

    gene_parsing_function = get_gene_parsing_function(germline_gene)

    gene_path = (
        Path(db_dir) / "gene_db" / f"{germline_gene.value.lower()}_genes" / organism.value / f"{locus.value}.fasta"
    )
    genes = read_genes(gene_path, gene_parsing_function)

    gene_params = params.nt_gene(germline_gene, locus)
    genes_aligner = GeneAligner(genes, **gene_params.kwargs(), ssw_params=gene_params.aligner)

    return genes_aligner


def create_aa_v_gene_aligner(
    allowed_species: Sequence[Organism], aa_genes_dir: Path, params: AlignerParams = DEFAULT_ALIGNER_PARAMS
) -> GeneAlignerAA:
    genes = []

    if not allowed_species:
        allowed_species = (Organism.HOMO_SAPIENS, Organism.MUS_MUSCULUS, Organism.VICUGNA_PACOS)

    for species in allowed_species:
        input_path = aa_genes_dir / "v_genes" / f"{species.value}.fasta"

        assert input_path.exists(), f"Input germline file: {str(input_path)} does not exists"

        species_genes = read_genes(input_path, partial(parse_gene_aa, species=species))

        genes.extend(species_genes)

    genes_aligner = GeneAlignerAA(genes, **params.v_aa.kwargs(), ssw_params=params.v_aa.aligner)

    return genes_aligner


def create_aa_j_gene_aligner(
    organism: Organism, locus: Locus, aa_genes_dir: Path, params: AlignerParams = DEFAULT_ALIGNER_PARAMS
) -> GeneAlignerAA:
    input_path = aa_genes_dir / "j_genes" / organism.value / f"{locus.value}.fasta"
    assert input_path.exists(), f"Input germline file: {str(input_path)} does not exists"

    genes = read_genes(input_path, partial(parse_gene_aa, species=organism))

    genes_aligner = GeneAlignerAA(genes, **params.j_aa.kwargs(), ssw_params=params.j_aa.aligner)

    return genes_aligner


def create_aa_c_gene_aligner(
    organism: Organism, locus: Locus, aa_genes_dir: Path, params: AlignerParams = DEFAULT_ALIGNER_PARAMS
) -> GeneAlignerAA:
    input_path = aa_genes_dir / "c_genes" / organism.value / f"{locus.value}.fasta"
    assert input_path.exists(), f"Input germline file: {str(input_path)} does not exists"

    genes = read_genes(input_path, partial(parse_gene_aa, species=organism))

    genes_aligner = GeneAlignerAA(genes, **params.c_aa.kwargs(), ssw_params=params.c_aa.aligner)

    return genes_aligner
