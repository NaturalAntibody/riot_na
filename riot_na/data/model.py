import json
from dataclasses import dataclass
from enum import Enum, unique
from typing import Generic, Optional, Sequence, TypeVar

from riot_na.common.serialization_utils import base64_encode


class Cigar(str):
    def __new__(cls, *args, **kwargs):
        return str.__new__(cls, *args, **kwargs)


class AlignmentString(str):
    def __new__(cls, *args, **kwargs):
        return str.__new__(cls, *args, **kwargs)


@dataclass(frozen=True)
class InternalAlignmentEntry:  # pylint: disable=too-many-instance-attributes
    """https://github.com/soedinglab/MMseqs2/wiki#internal-alignment-format"""

    target_id: str
    alignment_score: float
    seq_identity: float
    e_value: float
    q_start: int
    q_end: int
    q_len: int
    t_start: int
    t_end: int
    t_len: int
    cigar: Cigar
    query: str
    rev_comp: bool

    def __lt__(self, other):
        if self.e_value == other.e_value:
            return (self.t_end - self.t_start) < (other.t_end - other.t_start)
        return self.e_value < other.e_value

    def __gt__(self, other):
        if self.e_value == other.e_value:
            return (self.t_end - self.t_start) > (other.t_end - other.t_start)
        return self.e_value > other.e_value


@dataclass(frozen=True)
class InternalAlignmentEntryAA:  # pylint: disable=too-many-instance-attributes
    target_id: str
    alignment_score: Optional[float]
    seq_identity: Optional[float]
    e_value: Optional[float]
    q_start: int
    q_end: int
    t_start: int
    t_end: int
    cigar: Cigar

    def __lt__(self, other):
        if self.e_value == other.e_value:
            return (self.t_end - self.t_start) < (other.t_end - other.t_start)
        return self.e_value < other.e_value

    def __gt__(self, other):
        if self.e_value == other.e_value:
            return (self.t_end - self.t_start) > (other.t_end - other.t_start)
        return self.e_value > other.e_value


@dataclass(frozen=True)
class SchemeAlignment:
    target_id: str
    q_start: int
    q_end: int
    alignment_str: AlignmentString
    exc: Optional[str] = None


class Locus(str, Enum):
    IGH = "igh"
    IGL = "igl"
    IGK = "igk"


class ChainType(str, Enum):
    HEAVY = "H"
    LIGHT = "L"

    @classmethod
    def from_locus(cls, locus: Locus):
        match locus:
            case Locus.IGH:
                return cls.HEAVY
            case _:
                return cls.LIGHT


@unique
class Scheme(str, Enum):
    KABAT = "kabat"
    CHOTHIA = "chothia"
    IMGT = "imgt"
    MARTIN = "martin"


@unique
class Organism(str, Enum):
    HOMO_SAPIENS = "human"
    MUS_MUSCULUS = "mouse"
    VICUGNA_PACOS = "alpaca"
    CUSTOM = "custom"


class ShortRegion(str, Enum):
    FW1 = "fw1"
    FW2 = "fw2"
    FW3 = "fw3"
    FW4 = "fw4"
    CDR1 = "cdr1"
    CDR2 = "cdr2"
    CDR3 = "cdr3"


EntryValidationResults = dict[str, Optional[bool]]


@dataclass
class AirrRearrangementEntryNT:  # pylint: disable=too-many-instance-attributes
    """https://docs.airr-community.org/en/stable/datarep/rearrangements.html"""

    sequence_header: str
    sequence: str
    numbering_scheme: str
    locus: Optional[str] = None
    locus_species: Optional[str] = None
    stop_codon: Optional[bool] = None
    vj_in_frame: Optional[bool] = None
    v_frameshift: Optional[bool] = None
    j_frameshift: Optional[bool] = None
    productive: Optional[bool] = None
    rev_comp: Optional[bool] = None
    complete_vdj: Optional[bool] = None
    v_call: Optional[str] = None
    d_call: Optional[str] = None
    j_call: Optional[str] = None
    c_call: Optional[str] = None
    v_frame: Optional[int] = None
    j_frame: Optional[int] = None
    sequence_alignment: Optional[str] = None
    germline_alignment: Optional[str] = None
    sequence_aa: Optional[str] = None
    sequence_alignment_aa: Optional[str] = None
    germline_alignment_aa: Optional[str] = None
    v_alignment_start: Optional[int] = None
    v_alignment_end: Optional[int] = None
    d_alignment_start: Optional[int] = None
    d_alignment_end: Optional[int] = None
    j_alignment_start: Optional[int] = None
    j_alignment_end: Optional[int] = None
    v_sequence_alignment: Optional[str] = None
    v_sequence_alignment_aa: Optional[str] = None
    v_germline_alignment: Optional[str] = None
    v_germline_alignment_aa: Optional[str] = None
    d_sequence_alignment: Optional[str] = None

    d_germline_alignment: Optional[str] = None
    j_sequence_alignment: Optional[str] = None
    j_sequence_alignment_aa: Optional[str] = None
    j_germline_alignment: Optional[str] = None
    j_germline_alignment_aa: Optional[str] = None
    c_sequence_alignment: Optional[str] = None
    c_sequence_alignment_aa: Optional[str] = None
    c_germline_alignment: Optional[str] = None
    c_germline_alignment_aa: Optional[str] = None
    fwr1: Optional[str] = None
    fwr1_aa: Optional[str] = None
    cdr1: Optional[str] = None
    cdr1_aa: Optional[str] = None
    fwr2: Optional[str] = None
    fwr2_aa: Optional[str] = None
    cdr2: Optional[str] = None
    cdr2_aa: Optional[str] = None
    fwr3: Optional[str] = None
    fwr3_aa: Optional[str] = None
    cdr3: Optional[str] = None
    cdr3_aa: Optional[str] = None
    fwr4: Optional[str] = None
    fwr4_aa: Optional[str] = None
    junction: Optional[str] = None
    junction_aa: Optional[str] = None
    junction_length: Optional[int] = None
    junction_aa_length: Optional[int] = None
    v_score: Optional[float] = None
    d_score: Optional[float] = None
    j_score: Optional[float] = None
    c_score: Optional[float] = None

    v_cigar: Optional[str] = None
    d_cigar: Optional[str] = None
    j_cigar: Optional[str] = None
    c_cigar: Optional[str] = None

    v_support: Optional[float] = None
    d_support: Optional[float] = None
    j_support: Optional[float] = None
    c_support: Optional[float] = None
    v_identity: Optional[float] = None
    d_identity: Optional[float] = None
    j_identity: Optional[float] = None
    c_identity: Optional[float] = None
    v_sequence_start: Optional[int] = None
    v_sequence_end: Optional[int] = None
    d_sequence_start: Optional[int] = None
    d_sequence_end: Optional[int] = None
    j_sequence_start: Optional[int] = None
    j_sequence_end: Optional[int] = None
    c_sequence_start: Optional[int] = None
    c_sequence_end: Optional[int] = None

    v_germline_start: Optional[int] = None
    v_germline_end: Optional[int] = None
    d_germline_start: Optional[int] = None
    d_germline_end: Optional[int] = None
    j_germline_start: Optional[int] = None
    j_germline_end: Optional[int] = None
    c_germline_start: Optional[int] = None
    c_germline_end: Optional[int] = None

    fwr1_start: Optional[int] = None
    fwr1_end: Optional[int] = None
    cdr1_start: Optional[int] = None
    cdr1_end: Optional[int] = None
    fwr2_start: Optional[int] = None
    fwr2_end: Optional[int] = None
    cdr2_start: Optional[int] = None
    cdr2_end: Optional[int] = None
    fwr3_start: Optional[int] = None
    fwr3_end: Optional[int] = None
    cdr3_start: Optional[int] = None
    cdr3_end: Optional[int] = None
    fwr4_start: Optional[int] = None
    fwr4_end: Optional[int] = None

    sequence_aa_scheme_cigar: Optional[Cigar] = None
    scheme_residue_mapping: Optional[dict[str, str]] = None
    positional_scheme_mapping: Optional[dict[int, str]] = None
    exc: Optional[str] = None
    additional_validation_flags: Optional[EntryValidationResults] = None

    def __gt__(self, other):
        if self.v_score is None:
            return False
        if other.v_score is None:
            return True
        return self.v_score > other.v_score

    @property
    def primary_sequence_aa(self) -> str:
        if self.scheme_residue_mapping is None:
            return ""
        return "".join(self.scheme_residue_mapping.values())


@dataclass
class AirrRearrangementEntryAA:  # pylint: disable=too-many-instance-attributes
    """https://docs.airr-community.org/en/stable/datarep/rearrangements.html"""

    sequence_header: str
    sequence_aa: str
    numbering_scheme: str
    locus: Optional[str] = None
    locus_species: Optional[str] = None
    stop_codon: Optional[bool] = None
    productive: Optional[bool] = None  # V and J aligned and no stop codon
    complete_vdj: Optional[bool] = None
    v_call: Optional[str] = None
    j_call: Optional[str] = None
    c_call: Optional[str] = None

    germline_alignment_aa: Optional[str] = None
    sequence_alignment_aa: Optional[str] = None
    v_alignment_start_aa: Optional[int] = None
    v_alignment_end_aa: Optional[int] = None
    j_alignment_start_aa: Optional[int] = None
    j_alignment_end_aa: Optional[int] = None
    v_sequence_alignment_aa: Optional[str] = None
    v_germline_alignment_aa: Optional[str] = None
    j_sequence_alignment_aa: Optional[str] = None
    j_germline_alignment_aa: Optional[str] = None
    c_sequence_alignment_aa: Optional[str] = None
    c_germline_alignment_aa: Optional[str] = None

    fwr1_aa: Optional[str] = None
    cdr1_aa: Optional[str] = None
    fwr2_aa: Optional[str] = None
    cdr2_aa: Optional[str] = None
    fwr3_aa: Optional[str] = None
    cdr3_aa: Optional[str] = None
    fwr4_aa: Optional[str] = None
    junction_aa: Optional[str] = None
    junction_aa_length: Optional[int] = None

    v_score_aa: Optional[float] = None
    j_score_aa: Optional[float] = None
    c_score_aa: Optional[float] = None

    v_cigar_aa: Optional[str] = None
    j_cigar_aa: Optional[str] = None
    c_cigar_aa: Optional[str] = None

    v_support_aa: Optional[float] = None
    j_support_aa: Optional[float] = None
    c_support_aa: Optional[float] = None
    v_identity_aa: Optional[float] = None
    j_identity_aa: Optional[float] = None
    c_identity_aa: Optional[float] = None

    v_sequence_start_aa: Optional[int] = None
    v_sequence_end_aa: Optional[int] = None
    j_sequence_start_aa: Optional[int] = None
    j_sequence_end_aa: Optional[int] = None
    c_sequence_start_aa: Optional[int] = None
    c_sequence_end_aa: Optional[int] = None

    v_germline_start_aa: Optional[int] = None
    v_germline_end_aa: Optional[int] = None
    j_germline_start_aa: Optional[int] = None
    j_germline_end_aa: Optional[int] = None
    c_germline_start_aa: Optional[int] = None
    c_germline_end_aa: Optional[int] = None

    fwr1_start_aa: Optional[int] = None
    fwr1_end_aa: Optional[int] = None
    cdr1_start_aa: Optional[int] = None
    cdr1_end_aa: Optional[int] = None
    fwr2_start_aa: Optional[int] = None
    fwr2_end_aa: Optional[int] = None
    cdr2_start_aa: Optional[int] = None
    cdr2_end_aa: Optional[int] = None
    fwr3_start_aa: Optional[int] = None
    fwr3_end_aa: Optional[int] = None
    cdr3_start_aa: Optional[int] = None
    cdr3_end_aa: Optional[int] = None
    fwr4_start_aa: Optional[int] = None
    fwr4_end_aa: Optional[int] = None

    sequence_aa_scheme_cigar: Optional[Cigar] = None
    scheme_residue_mapping: Optional[dict[str, str]] = None
    positional_scheme_mapping: Optional[dict[int, str]] = None
    exc: Optional[str] = None
    additional_validation_flags: Optional[EntryValidationResults] = None

    def __gt__(self, other):
        if self.v_score_aa is None:
            return False
        if other.v_score_aa is None:
            return True
        return self.v_score_aa > other.v_score_aa

    @property
    def primary_sequence_aa(self) -> str:
        if self.scheme_residue_mapping is None:
            return ""
        return "".join(self.scheme_residue_mapping.values())


@dataclass
class SegmentedAirrRearrangementEntryAA(AirrRearrangementEntryAA):
    query_sequence: Optional[str] = None
    segment_start: Optional[int] = None
    segment_end: Optional[int] = None


@dataclass
class SegmentedAirrRearrangementEntryNT(AirrRearrangementEntryNT):
    query_sequence: Optional[str] = None
    segment_start: Optional[int] = None
    segment_end: Optional[int] = None


AirrRearrangementEntry_co = TypeVar(
    "AirrRearrangementEntry_co",
    AirrRearrangementEntryNT,
    AirrRearrangementEntryAA,
    SegmentedAirrRearrangementEntryNT,
    SegmentedAirrRearrangementEntryAA,
    covariant=True,
)


def serialize_airr_entry(entry: AirrRearrangementEntry_co) -> dict[str, str]:
    entry_dict = entry.__dict__
    entry_dict["scheme_residue_mapping"] = entry_dict["scheme_residue_mapping"] and base64_encode(
        json.dumps(entry_dict["scheme_residue_mapping"], separators=(",", ":"))
    )
    entry_dict["positional_scheme_mapping"] = entry_dict["positional_scheme_mapping"] and base64_encode(
        json.dumps(entry_dict["positional_scheme_mapping"], separators=(",", ":"))
    )
    entry_dict["exc"] = entry_dict["exc"] and base64_encode(entry_dict["exc"])
    entry_dict["additional_validation_flags"] = entry_dict["additional_validation_flags"] and base64_encode(
        json.dumps(entry_dict["additional_validation_flags"], separators=(",", ":"))
    )

    return entry_dict


@dataclass
class RegionOffsetsAA:  # pylint: disable=too-many-instance-attributes
    fwr1_start_aa: int
    fwr1_end_aa: int
    cdr1_start_aa: int
    cdr1_end_aa: int
    fwr2_start_aa: int
    fwr2_end_aa: int
    cdr2_start_aa: int
    cdr2_end_aa: int
    fwr3_start_aa: int
    fwr3_end_aa: int
    cdr3_start_aa: int
    cdr3_end_aa: int
    fwr4_start_aa: int
    fwr4_end_aa: int


@dataclass
class RegionOffsetsNT:  # pylint: disable=too-many-instance-attributes
    fwr1_start: int
    fwr1_end: int
    cdr1_start: int
    cdr1_end: int
    fwr2_start: int
    fwr2_end: int
    cdr2_start: int
    cdr2_end: int
    fwr3_start: int
    fwr3_end: int
    cdr3_start: int
    cdr3_end: int
    fwr4_start: int
    fwr4_end: int


GeneId = str
GeneSequence = str


@dataclass(frozen=True)
class Gene:
    name: GeneId
    species: Organism
    locus: Locus
    sequence: GeneSequence
    reading_frame: Optional[int]


@dataclass(frozen=True)
class GeneAA:
    name: GeneId
    species: Organism
    locus: Locus
    sequence: GeneSequence


@dataclass(frozen=True)
class SpeciesGeneMatch:
    species_gene_id: str
    gene_id: str
    rev_comp: bool
    coverage: int
    species: Organism
    locus: Locus


@dataclass(frozen=True)
class SpeciesPrefilteringResult:
    query: str
    rev_comp_query: str
    top_matches: list[SpeciesGeneMatch]


@dataclass(frozen=True)
class SegmentMatch:
    segment_start: int
    segment_length: int
    query_start: int
    query_end: int
    coverage: int
    match_count: int
    matching_genes: list[SpeciesGeneMatch]


@dataclass(frozen=True)
class SpeciesPrefilteringSegmentResult:
    query: str
    rev_comp_query: str
    segments: list[SegmentMatch]


@dataclass(frozen=True)
class AlignmentEntryNT:  # pylint: disable=too-many-instance-attributes
    target_id: str
    alignment_score: float
    seq_identity: float
    e_value: float
    q_start: int
    q_end: int
    q_len: int
    t_start: int
    t_end: int
    t_len: int
    cigar: Cigar
    rev_comp: bool

    species: Organism
    locus: Locus
    q_seq: GeneSequence
    t_seq: GeneSequence
    reading_frame: Optional[int]

    def __lt__(self, other):
        if self.e_value == other.e_value:
            if self.e_value == 0:
                if self.seq_identity == other.seq_identity:
                    return self.target_id < other.target_id
                return self.seq_identity > other.seq_identity
            return self.target_id < other.target_id
        return self.e_value < other.e_value

    def __gt__(self, other):
        if self.e_value == other.e_value:
            if self.e_value == 0:
                if self.seq_identity == other.seq_identity:
                    return self.target_id > other.target_id
                return self.seq_identity < other.seq_identity
            return self.target_id > other.target_id
        return self.e_value > other.e_value

    def lookup_gene_id(self) -> str:
        return f"{self.species.value}|{self.locus.value}|{self.target_id}"


@dataclass(frozen=True)
class AlignmentEntryAA:  # pylint: disable=too-many-instance-attributes
    target_id: str
    alignment_score: Optional[float]
    seq_identity: Optional[float]
    e_value: Optional[float]
    q_start: int
    q_end: int
    t_start: int
    t_end: int
    cigar: Cigar

    species: Organism
    locus: Locus
    q_seq: GeneSequence
    t_seq: GeneSequence

    def __lt__(self, other):
        if self.e_value == other.e_value:
            if self.e_value == 0:
                if self.seq_identity == other.seq_identity:
                    return self.target_id < other.target_id
                return self.seq_identity > other.seq_identity
            return self.target_id < other.target_id
        return self.e_value < other.e_value

    def __gt__(self, other):
        if self.e_value == other.e_value:
            if self.e_value == 0:
                if self.seq_identity == other.seq_identity:
                    return self.target_id > other.target_id
                return self.seq_identity < other.seq_identity
            return self.target_id > other.target_id
        return self.e_value > other.e_value

    def lookup_gene_id(self) -> str:
        return f"{self.species.value}|{self.locus.value}|{self.target_id}"


AlignmentEntry = TypeVar("AlignmentEntry", AlignmentEntryNT, AlignmentEntryAA)


@dataclass
class AlignmentSegment(Generic[AlignmentEntry]):
    alignments: Sequence[AlignmentEntry]
    best_alignment: Optional[AlignmentEntry] = None  # Best alignment for this segment
    segment_start: Optional[int] = None
    segment_end: Optional[int] = None

    def __post_init__(self):
        """Set the best alignment after initialization"""
        if self.alignments:
            self.best_alignment = min(self.alignments)  # Lowest e-value first


@dataclass
class AlignmentsNT:
    v: Optional[AlignmentEntryNT] = None
    d: Optional[AlignmentEntryNT] = None
    j: Optional[AlignmentEntryNT] = None
    c: Optional[AlignmentEntryNT] = None
    segment_start: Optional[int] = None
    segment_end: Optional[int] = None


@dataclass
class AlignmentsAA:
    v: Optional[AlignmentEntryAA] = None
    j: Optional[AlignmentEntryAA] = None
    c: Optional[AlignmentEntryAA] = None
    segment_start: Optional[int] = None
    segment_end: Optional[int] = None


@dataclass
class TranslatedAlignmentsNT:
    translated_query: str
    reading_frame: int
    aa_alignments: AlignmentsAA


class GermlineGene(Enum):
    V = "V"
    D = "D"
    J = "J"
    C = "C"


class InputType(str, Enum):
    NT = "nt"
    AA = "aa"
