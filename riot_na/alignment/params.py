"""
Thresholds and Smith-Waterman parameters of the gene aligners.

All dataclasses are frozen (hence hashable, so they can be used as arguments of the cached ``get_or_create_*``
factories). To change a setting create a modified copy and pass it as the last argument of a factory::

    from dataclasses import replace
    from riot_na import DEFAULT_ALIGNER_PARAMS, create_riot_nt

    params = replace(DEFAULT_ALIGNER_PARAMS, ch_nt=replace(DEFAULT_ALIGNER_PARAMS.ch_nt, e_value_threshold=1e-3))
    riot = create_riot_nt(params=params)
"""

from __future__ import annotations

import warnings
from dataclasses import asdict, dataclass, replace
from typing import Any, Optional

from riot_na.data.constants import BLOSUM_62
from riot_na.data.model import GermlineGene, Locus, Organism


@dataclass(frozen=True, kw_only=True)
class NTSSWParams:
    """Striped Smith-Waterman parameters of nucleotide alignments. A gap of length n costs open + (n - 1) * extend."""

    match_score: int = 1
    mismatch_score: int = -1
    gap_open_penalty: int = 4
    gap_extend_penalty: int = 1

    def kwargs(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, kw_only=True)
class AASSWParams:
    """Striped Smith-Waterman parameters (BLOSUM62) of amino acid alignments."""

    gap_open_penalty: int = 11
    gap_extend_penalty: int = 1

    def kwargs(self) -> dict[str, Any]:
        return {**asdict(self), "protein": True, "substitution_matrix": BLOSUM_62}


@dataclass(frozen=True, kw_only=True)
class GeneParams:
    """
    Prefiltering and filtering thresholds of a gene aligner.

    - kmer_size, distance_threshold, top_n, modulo_n: prefiltering of germlines (see MultiSpeciesPrefiltering).
    - e_value_threshold: alignments with a higher E-value are dropped.
    - alignment_length_threshold: alignments with a shorter germline part are dropped.
    - min_seq_identity: alignments with a lower sequence identity (0-1) are dropped.
    """

    kmer_size: int
    distance_threshold: int
    top_n: int
    modulo_n: int
    e_value_threshold: Optional[float] = None
    alignment_length_threshold: Optional[int] = None
    min_seq_identity: Optional[float] = None

    def kwargs(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, kw_only=True)
class AAGeneParams(GeneParams):
    """
    Amino acid gene aligner thresholds.

    - max_cdr3_length: longest allowed CDR3 (J alignment).
    - min_internal_alignment_length: alignments shorter than that are dropped, unless they are fragments at the end
      of the germline.
    """

    max_cdr3_length: Optional[int] = None
    min_internal_alignment_length: Optional[int] = None


@dataclass(frozen=True, kw_only=True)
class SegmentGeneParams(GeneParams):
    """Thresholds of the V aligner finding domains of multiple domain nucleotide sequences."""

    min_prefiltering_coverage: int
    min_segment_length: int


@dataclass(frozen=True, kw_only=True)
class SegmentAAGeneParams(SegmentGeneParams):
    """Thresholds of the V aligner finding domains of multiple domain amino acid sequences."""

    max_cdr3_length: Optional[int] = None


@dataclass(frozen=True, kw_only=True)
class AlignerParams:  # pylint: disable=too-many-instance-attributes
    """
    Thresholds and Smith-Waterman parameters of all gene aligners. The defaults are the established RIOT settings.

    Nucleotide aligners: V, D, J by locus (jh, jk, jl) and C by locus (ch, ck, cl). Amino acid aligners: V, J, C.
    ``segment_v_*`` are the V aligners finding domains of multiple domain sequences, the other genes of those are
    aligned with the same settings as of single domain sequences.
    """

    v_nt: GeneParams = GeneParams(kmer_size=5, distance_threshold=13, top_n=12, modulo_n=2, e_value_threshold=5.9e-08, alignment_length_threshold=100)  # fmt: skip
    d_nt: GeneParams = GeneParams(kmer_size=5, distance_threshold=3, top_n=5, modulo_n=1)
    jh_nt: GeneParams = GeneParams(kmer_size=5, distance_threshold=3, top_n=5, modulo_n=2)
    jk_nt: GeneParams = GeneParams(kmer_size=5, distance_threshold=5, top_n=5, modulo_n=1)
    jl_nt: GeneParams = GeneParams(kmer_size=5, distance_threshold=7, top_n=5, modulo_n=2)
    ch_nt: GeneParams = GeneParams(kmer_size=5, distance_threshold=5, top_n=5, modulo_n=2, e_value_threshold=1e-05, min_seq_identity=0.85)  # fmt: skip
    ck_nt: GeneParams = GeneParams(kmer_size=3, distance_threshold=5, top_n=3, modulo_n=1, e_value_threshold=1e-05, min_seq_identity=0.85)  # fmt: skip
    cl_nt: GeneParams = GeneParams(kmer_size=3, distance_threshold=5, top_n=3, modulo_n=1, e_value_threshold=1e-05, min_seq_identity=0.85)  # fmt: skip

    v_aa: AAGeneParams = AAGeneParams(kmer_size=3, distance_threshold=3, top_n=20, modulo_n=1, e_value_threshold=6.79e-15, alignment_length_threshold=50)  # fmt: skip
    j_aa: AAGeneParams = AAGeneParams(kmer_size=3, distance_threshold=3, top_n=5, modulo_n=1, e_value_threshold=0.05, alignment_length_threshold=5, max_cdr3_length=60)  # fmt: skip
    # C alignments of at least min_internal_alignment_length residues have E-values of at most about 1e-17
    # (therapeutics), with no candidates between that and 1. The E-value threshold lies in that gap, with margin for
    # C regions more mutated than in therapeutics. The minimum sequence identity filters out matches of non antibody
    # regions following the V region, e.g. proline rich CAR hinges to IgG and IgA hinges.
    c_aa: AAGeneParams = AAGeneParams(kmer_size=3, distance_threshold=3, top_n=5, modulo_n=2, e_value_threshold=1e-10, alignment_length_threshold=4, min_seq_identity=0.85, min_internal_alignment_length=50)  # fmt: skip

    segment_v_nt: SegmentGeneParams = SegmentGeneParams(kmer_size=5, distance_threshold=13, top_n=12, modulo_n=2, e_value_threshold=5.8e-08, alignment_length_threshold=100, min_prefiltering_coverage=75, min_segment_length=100)  # fmt: skip
    segment_v_aa: SegmentAAGeneParams = SegmentAAGeneParams(kmer_size=3, distance_threshold=4, top_n=12, modulo_n=1, e_value_threshold=1.0313e-18, alignment_length_threshold=50, min_prefiltering_coverage=32, min_segment_length=35)  # fmt: skip

    nt_ssw: NTSSWParams = NTSSWParams()
    aa_ssw: AASSWParams = AASSWParams()
    # C regions of engineered antibodies can be separated from the preceding domains by long linkers. A high gap
    # extension penalty keeps the local alignment from bridging them, e.g. a (GGGGS)n linker between two Fc copies,
    # to align the preceding domain to the C germline.
    c_aa_ssw: AASSWParams = AASSWParams(gap_extend_penalty=5)

    def nt_gene(self, germline_gene: GermlineGene, locus: Optional[Locus]) -> GeneParams:
        """Thresholds of the nucleotide aligner of the gene (and locus, for J and C)."""
        match germline_gene:
            case GermlineGene.V:
                return self.v_nt
            case GermlineGene.D:
                return self.d_nt
            case GermlineGene.J:
                assert locus is not None
                return {Locus.IGH: self.jh_nt, Locus.IGK: self.jk_nt, Locus.IGL: self.jl_nt}[locus]
            case GermlineGene.C:
                assert locus is not None
                return {Locus.IGH: self.ch_nt, Locus.IGK: self.ck_nt, Locus.IGL: self.cl_nt}[locus]
        raise ValueError(f"Nucleotide aligner not supported for {germline_gene}")

    def aa_gene(self, germline_gene: GermlineGene) -> AAGeneParams:
        """Thresholds of the amino acid aligner of the gene."""
        match germline_gene:
            case GermlineGene.V:
                return self.v_aa
            case GermlineGene.J:
                return self.j_aa
            case GermlineGene.C:
                return self.c_aa
        raise ValueError("AA aligner not supported for gene type")

    def with_c_min_seq_identity(self, min_seq_identity: Optional[float]) -> AlignerParams:
        """Copy with the minimum sequence identity of constant region alignments (nucleotide and amino acid) set."""
        return replace(
            self,
            ch_nt=replace(self.ch_nt, min_seq_identity=min_seq_identity),
            ck_nt=replace(self.ck_nt, min_seq_identity=min_seq_identity),
            cl_nt=replace(self.cl_nt, min_seq_identity=min_seq_identity),
            c_aa=replace(self.c_aa, min_seq_identity=min_seq_identity),
        )


DEFAULT_ALIGNER_PARAMS = AlignerParams()


def warn_c_genes_organism(c_genes_organism: Optional[Organism], params: AlignerParams) -> None:
    """Warn that constant regions are aligned to germlines of another organism with the sequence identity filter."""
    if c_genes_organism is None:
        return
    nt_identity = params.ch_nt.min_seq_identity
    aa_identity = params.c_aa.min_seq_identity
    if nt_identity is None and aa_identity is None:
        return

    def percent(value: Optional[float]) -> str:
        return "no limit" if value is None else f"{value:.0%}"

    warnings.warn(
        f"c_genes_organism={c_genes_organism.name}: constant regions are aligned to {c_genes_organism.value} germlines "
        f"for V alignments of any species. The current minimum sequence identity of C alignments is "
        f"{percent(nt_identity)} (nucleotide) and {percent(aa_identity)} (amino acid). Constant regions of another "
        "organism are less similar to the germlines, so you should probably lower it. To do so create a modified "
        "params object and pass it as the last argument, e.g.: "
        "`params = DEFAULT_ALIGNER_PARAMS.with_c_min_seq_identity(0.4); create_riot_nt(c_genes_organism=..., "
        "params=params)` (see riot_na.alignment.params to change other thresholds).",
        UserWarning,
        stacklevel=3,
    )
