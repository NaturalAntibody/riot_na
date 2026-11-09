"""
Thresholds and Smith-Waterman parameters of the gene aligners.

All dataclasses are frozen (hence hashable, so they can be used as arguments of the cached ``get_or_create_*``
factories). To change a setting create a modified copy and pass it as the last argument of a factory::

    from dataclasses import replace
    from riot_na import DEFAULT_ALIGNER_PARAMS, create_riot_nt

    params = DEFAULT_ALIGNER_PARAMS.with_param("ch_nt.thresholds.e_value_threshold", 1e-3)
    riot = create_riot_nt(params=params)

Every gene aligner has ``prefiltering`` (germline selection), ``aligner`` (Smith-Waterman scoring) and ``thresholds``
(what alignments must pass) settings.
"""

from __future__ import annotations

import warnings
from dataclasses import asdict, dataclass, fields, is_dataclass, replace
from typing import Any, Optional, Union

from riot_na.data.constants import BLOSUM_62, BLOSUM_80
from riot_na.data.model import GermlineGene, Locus, Organism

SUBSTITUTION_MATRICES: dict[str, dict[str, dict[str, int]]] = {"blosum62": BLOSUM_62, "blosum80": BLOSUM_80}


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
    """
    Striped Smith-Waterman parameters of amino acid alignments. ``matrix`` is the name of the substitution matrix
    (see SUBSTITUTION_MATRICES). A gap of length n costs open + (n - 1) * extend.
    """

    gap_open_penalty: int = 11
    gap_extend_penalty: int = 1
    matrix: str = "blosum62"

    def substitution_matrix(self) -> dict[str, dict[str, int]]:
        try:
            return SUBSTITUTION_MATRICES[self.matrix]
        except KeyError:
            raise ValueError(
                f"Unknown substitution matrix {self.matrix!r}, use one of {list(SUBSTITUTION_MATRICES)}"
            ) from None

    def kwargs(self) -> dict[str, Any]:
        return {
            "gap_open_penalty": self.gap_open_penalty,
            "gap_extend_penalty": self.gap_extend_penalty,
            "protein": True,
            "substitution_matrix": self.substitution_matrix(),
        }


@dataclass(frozen=True, kw_only=True)
class PrefilteringParams:
    """
    Prefiltering of germlines by k-mers, which picks the germlines to align to (see MultiSpeciesPrefiltering).

    - kmer_size: length of the k-mers.
    - distance_threshold, top_n, modulo_n: selection of the best matching germlines.
    """

    kmer_size: int
    distance_threshold: int
    top_n: int
    modulo_n: int


@dataclass(frozen=True, kw_only=True)
class SegmentPrefilteringParams(PrefilteringParams):
    """Prefiltering of the V aligner finding domains of multiple domain sequences."""

    min_prefiltering_coverage: int
    min_segment_length: int


@dataclass(frozen=True, kw_only=True)
class Thresholds:
    """
    Thresholds the alignments to the prefiltered germlines must pass.

    - e_value_threshold: alignments with a higher E-value are dropped.
    - alignment_length_threshold: alignments with a shorter germline part are dropped.
    - min_seq_identity: alignments with a lower sequence identity (0-1) are dropped.
    """

    e_value_threshold: Optional[float] = None
    alignment_length_threshold: Optional[int] = None
    min_seq_identity: Optional[float] = None


@dataclass(frozen=True, kw_only=True)
class AAThresholds(Thresholds):
    """
    Filters of amino acid alignments.

    - max_cdr3_length: longest allowed CDR3 (J alignment).
    - min_internal_alignment_length: alignments shorter than that are dropped, unless they are fragments at the end
      of the germline.
    """

    max_cdr3_length: Optional[int] = None
    min_internal_alignment_length: Optional[int] = None


@dataclass(frozen=True, kw_only=True)
class SegmentAAThresholds(Thresholds):
    """Filters of amino acid alignments of the V aligner finding domains of multiple domain sequences."""

    max_cdr3_length: Optional[int] = None


@dataclass(frozen=True, kw_only=True)
class GeneParams:
    """
    Settings of a nucleotide gene aligner.

    - prefiltering: selection of the germlines to align to.
    - aligner: Smith-Waterman scoring parameters.
    - thresholds: what an alignment must pass to become a candidate.
    """

    prefiltering: PrefilteringParams
    aligner: NTSSWParams = NTSSWParams()
    thresholds: Thresholds = Thresholds()

    def kwargs(self) -> dict[str, Any]:
        """Prefiltering and threshold settings as the (flat) keyword arguments of the gene aligner."""
        return {**asdict(self.prefiltering), **asdict(self.thresholds)}


@dataclass(frozen=True, kw_only=True)
class AAGeneParams(GeneParams):
    """Settings of an amino acid gene aligner."""

    aligner: AASSWParams = AASSWParams()  # type: ignore[assignment]
    thresholds: AAThresholds = AAThresholds()  # type: ignore[assignment]


@dataclass(frozen=True, kw_only=True)
class SegmentGeneParams(GeneParams):
    """Settings of the V aligner finding domains of multiple domain nucleotide sequences."""

    prefiltering: SegmentPrefilteringParams  # type: ignore[assignment]


@dataclass(frozen=True, kw_only=True)
class SegmentAAGeneParams(GeneParams):
    """Settings of the V aligner finding domains of multiple domain amino acid sequences."""

    prefiltering: SegmentPrefilteringParams  # type: ignore[assignment]
    aligner: AASSWParams = AASSWParams()  # type: ignore[assignment]
    thresholds: SegmentAAThresholds = SegmentAAThresholds()  # type: ignore[assignment]


@dataclass(frozen=True, kw_only=True)
class AlignerParams:  # pylint: disable=too-many-instance-attributes
    """
    Settings of all gene aligners (prefiltering, Smith-Waterman parameters and thresholds). The defaults are the established RIOT settings.

    Nucleotide aligners: V, D, J by locus (jh, jk, jl) and C by locus (ch, ck, cl). Amino acid aligners: V, J, C.
    ``segment_v_*`` are the V aligners finding domains of multiple domain sequences, the other genes of those are
    aligned with the same settings as of single domain sequences.
    """

    v_nt: GeneParams = GeneParams(
        prefiltering=PrefilteringParams(kmer_size=5, distance_threshold=13, top_n=12, modulo_n=2),
        thresholds=Thresholds(e_value_threshold=9.0e-08, alignment_length_threshold=100),
    )
    d_nt: GeneParams = GeneParams(
        prefiltering=PrefilteringParams(kmer_size=5, distance_threshold=3, top_n=5, modulo_n=1)
    )
    jh_nt: GeneParams = GeneParams(
        prefiltering=PrefilteringParams(kmer_size=5, distance_threshold=3, top_n=5, modulo_n=2)
    )
    jk_nt: GeneParams = GeneParams(
        prefiltering=PrefilteringParams(kmer_size=5, distance_threshold=5, top_n=5, modulo_n=1)
    )
    jl_nt: GeneParams = GeneParams(
        prefiltering=PrefilteringParams(kmer_size=5, distance_threshold=7, top_n=5, modulo_n=2)
    )
    ch_nt: GeneParams = GeneParams(
        prefiltering=PrefilteringParams(kmer_size=5, distance_threshold=5, top_n=5, modulo_n=2),
        thresholds=Thresholds(e_value_threshold=1.35e-05, min_seq_identity=0.8),
    )
    ck_nt: GeneParams = GeneParams(
        prefiltering=PrefilteringParams(kmer_size=3, distance_threshold=5, top_n=3, modulo_n=1),
        thresholds=Thresholds(e_value_threshold=1.35e-05, min_seq_identity=0.8),
    )
    cl_nt: GeneParams = GeneParams(
        prefiltering=PrefilteringParams(kmer_size=3, distance_threshold=5, top_n=3, modulo_n=1),
        thresholds=Thresholds(e_value_threshold=1.35e-05, min_seq_identity=0.8),
    )

    v_aa: AAGeneParams = AAGeneParams(
        prefiltering=PrefilteringParams(kmer_size=3, distance_threshold=3, top_n=20, modulo_n=1),
        thresholds=AAThresholds(e_value_threshold=1.2e-13, alignment_length_threshold=50),
    )
    j_aa: AAGeneParams = AAGeneParams(
        prefiltering=PrefilteringParams(kmer_size=3, distance_threshold=3, top_n=5, modulo_n=1),
        thresholds=AAThresholds(e_value_threshold=0.3, alignment_length_threshold=5, max_cdr3_length=60),
    )
    # C alignments of at least min_internal_alignment_length residues have E-values of at most about 1e-17
    # (therapeutics), with no candidates between that and 1. The E-value threshold lies in that gap, with margin for
    # C regions more mutated than in therapeutics. The minimum sequence identity filters out matches of non antibody
    # regions following the V region, e.g. proline rich CAR hinges to IgG and IgA hinges.
    # C regions of engineered antibodies can be separated from the preceding domains by long linkers. A high gap
    # extension penalty keeps the local alignment from bridging them, e.g. a (GGGGS)n linker between two Fc copies,
    # to align the preceding domain to the C germline.
    c_aa: AAGeneParams = AAGeneParams(
        prefiltering=PrefilteringParams(kmer_size=3, distance_threshold=3, top_n=5, modulo_n=2),
        aligner=AASSWParams(gap_extend_penalty=5),
        thresholds=AAThresholds(
            e_value_threshold=4e-12,
            alignment_length_threshold=4,
            min_seq_identity=0.8,
            min_internal_alignment_length=50,
        ),
    )

    segment_v_nt: SegmentGeneParams = SegmentGeneParams(
        prefiltering=SegmentPrefilteringParams(
            kmer_size=5,
            distance_threshold=13,
            top_n=12,
            modulo_n=2,
            min_prefiltering_coverage=75,
            min_segment_length=100,
        ),
        thresholds=Thresholds(e_value_threshold=9.0e-08, alignment_length_threshold=100),
    )
    segment_v_aa: SegmentAAGeneParams = SegmentAAGeneParams(
        prefiltering=SegmentPrefilteringParams(
            kmer_size=3, distance_threshold=4, top_n=12, modulo_n=1, min_prefiltering_coverage=32, min_segment_length=35
        ),
        thresholds=SegmentAAThresholds(e_value_threshold=3.7e-17, alignment_length_threshold=50),
    )

    def ssw_params(self) -> list[Union[NTSSWParams, AASSWParams]]:
        """The distinct Smith-Waterman parameters of all gene aligners."""
        found = [getattr(self, f.name).aligner for f in fields(self) if isinstance(getattr(self, f.name), GeneParams)]
        return list(dict.fromkeys(found))

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

    def with_param(self, path: str, value: Any) -> AlignerParams:
        """
        Copy with the setting at the dotted ``path`` replaced, e.g.
        ``DEFAULT_ALIGNER_PARAMS.with_param("v_aa.prefiltering.kmer_size", 4)``. Calls can be chained.
        """
        return _replace_path(self, path.split("."), value, path)

    def with_c_min_seq_identity(self, min_seq_identity: Optional[float]) -> AlignerParams:
        """Copy with the minimum sequence identity of constant region alignments (nucleotide and amino acid) set."""

        def with_identity(params: Any) -> Any:
            return replace(params, thresholds=replace(params.thresholds, min_seq_identity=min_seq_identity))

        return replace(
            self,
            ch_nt=with_identity(self.ch_nt),
            ck_nt=with_identity(self.ck_nt),
            cl_nt=with_identity(self.cl_nt),
            c_aa=with_identity(self.c_aa),
        )


def _replace_path(obj: Any, names: list[str], value: Any, path: str) -> Any:
    name, *rest = names
    if not is_dataclass(obj) or name not in {f.name for f in fields(obj)}:
        valid = [f.name for f in fields(obj)] if is_dataclass(obj) else []
        raise AttributeError(f"Invalid setting {path!r}: {type(obj).__name__} has no field {name!r}, valid: {valid}")
    new_value = _replace_path(getattr(obj, name), rest, value, path) if rest else value
    return replace(obj, **{name: new_value})  # type: ignore[type-var]


DEFAULT_ALIGNER_PARAMS = AlignerParams()


def warn_c_genes_organism(c_genes_organism: Optional[Organism], params: AlignerParams) -> None:
    """
    Warn that constant regions are aligned to germlines of another organism with the default sequence identity filter.
    No warning if the identity limit was lowered (or removed) for all C aligners.
    """
    if c_genes_organism is None:
        return
    nt_identity = params.ch_nt.thresholds.min_seq_identity
    aa_identity = params.c_aa.thresholds.min_seq_identity

    def lowered(gene: str) -> bool:
        current = getattr(params, gene).thresholds.min_seq_identity
        default = getattr(DEFAULT_ALIGNER_PARAMS, gene).thresholds.min_seq_identity
        return current is None or (default is not None and current < default)

    if all(lowered(gene) for gene in ("ch_nt", "ck_nt", "cl_nt", "c_aa")):
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
        "params=params)` (see riot_na.thresholds.params to change other thresholds).",
        UserWarning,
        stacklevel=3,
    )
