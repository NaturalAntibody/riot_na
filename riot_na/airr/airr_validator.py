from enum import Enum
from functools import cache, partial
from typing import Callable, Literal, Optional

from riot_na.alignment.alignment_utils import translate, unfold_cigar
from riot_na.data.model import (
    AirrRearrangementEntry_co,
    AirrRearrangementEntryAA,
    AirrRearrangementEntryNT,
    ChainType,
    EntryValidationResults,
    Locus,
    Scheme,
)
from riot_na.data.scheme_definitions import get_indels_positions

Region = Literal["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]


class Gene(Enum):
    V = "v_"
    J = "j_"
    C = "c_"
    _ = ""


REGIONS: list[Region] = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]
REGION_TO_GENE: dict[Region, dict[Literal["start", "end"], Gene]] = {
    "fwr1": {"start": Gene.V, "end": Gene.V},
    "cdr1": {"start": Gene.V, "end": Gene.V},
    "fwr2": {"start": Gene.V, "end": Gene.V},
    "cdr2": {"start": Gene.V, "end": Gene.V},
    "fwr3": {"start": Gene.V, "end": Gene.V},
    "cdr3": {"start": Gene.V, "end": Gene.J},
    "fwr4": {"start": Gene.J, "end": Gene.J},
}

GENE_TO_REGIONS_START_END: dict[Gene, list[str]] = {
    Gene.V: [
        "fwr1_start",
        "fwr1_end",
        "cdr1_start",
        "cdr1_end",
        "fwr2_start",
        "fwr2_end",
        "cdr2_start",
        "cdr2_end",
        "fwr3_start",
        "fwr3_end",
        "cdr3_start",
    ],
    Gene.J: ["cdr3_end", "fwr4_start", "fwr4_end"],
}

IMGT_CONSERVED_RESIDUES_HEAVY = {
    "23": "C",
    "41": "W",
    "104": "C",
    "118_heavy": "W",
}
IMGT_CONSERVED_RESIDUES_LIGHT = {
    "23": "C",
    "41": "W",
    "104": "C",
    "118_light": "F",
}


class SeqType(Enum):
    NT = ""
    AA = "_aa"


KABAT_LIGHT_FRAMEWORK_INSERTIONS = {"106.1"}
KABAT_HEAVY_FRAMEWORK_INSERTIONS = {"82.1", "82.2", "82.3"}

CHOTHIA_LIGHT_FRAMEWORK_INSERTIONS = {"106.1"}
CHOTHIA_HEAVY_FRAMEWORK_INSERTIONS = {"82.1", "82.2", "82.3"}

MARTIN_LIGHT_FRAMEWORK_INSERTIONS = {"40.1", "68.1", "68.2", "68.3", "68.4", "68.5", "68.6", "68.7", "68.8", "107.1"}
MARTIN_HEAVY_FRAMEWORK_INSERTIONS = {"8.1", "8.2", "8.3", "8.4", "72.1", "72.2", "72.3"}


@cache
def get_framework_insertion_positions(scheme: Scheme, chain_type: ChainType) -> set[str]:
    # this is scheme specific
    match (scheme, chain_type):
        case Scheme.IMGT, _:
            return set()
        case Scheme.KABAT, ChainType.HEAVY:
            return KABAT_HEAVY_FRAMEWORK_INSERTIONS
        case Scheme.KABAT, ChainType.LIGHT:
            return KABAT_LIGHT_FRAMEWORK_INSERTIONS
        case Scheme.CHOTHIA, ChainType.HEAVY:
            return CHOTHIA_HEAVY_FRAMEWORK_INSERTIONS
        case Scheme.CHOTHIA, ChainType.LIGHT:
            return CHOTHIA_LIGHT_FRAMEWORK_INSERTIONS
        case Scheme.MARTIN, ChainType.HEAVY:
            return MARTIN_HEAVY_FRAMEWORK_INSERTIONS
        case Scheme.MARTIN, ChainType.LIGHT:
            return MARTIN_LIGHT_FRAMEWORK_INSERTIONS

    raise ValueError(f"Unknown scheme {scheme}")


def _getattr(rearrangement: AirrRearrangementEntry_co, attr: str, seq_type: SeqType = SeqType.NT, gene: Gene = Gene._):
    return getattr(rearrangement, f"{gene.value}{attr}{seq_type.value}")


def validate_regions_in_aligned_sequence(
    rearrangement: AirrRearrangementEntry_co, seq_type: SeqType = SeqType.NT
) -> Optional[bool]:
    concatenated_regions = ""
    for region in REGIONS:
        if not _getattr(rearrangement, region, seq_type):
            continue
        concatenated_regions += _getattr(rearrangement, region, seq_type)

    sequence_alignment = _getattr(rearrangement, "sequence_alignment", seq_type).replace("-", "")
    if (not sequence_alignment) or (not concatenated_regions):
        return None
    return concatenated_regions in sequence_alignment


def validate_translated_regions_in_aligned_sequence_aa(rearrangement: AirrRearrangementEntryNT) -> Optional[bool]:
    concatenated_regions = ""
    for region in REGIONS:
        if not _getattr(rearrangement, region):
            continue
        concatenated_regions += _getattr(rearrangement, region)

    sequence_alignment = rearrangement.sequence_alignment_aa
    v_frame = rearrangement.v_frame
    if (not sequence_alignment) or (not concatenated_regions) or (v_frame is None):
        return None
    translated_regions = translate(concatenated_regions, 0)
    return translated_regions in sequence_alignment.replace("-", "")


def validate_consecutive_offsets(
    rearrangement: AirrRearrangementEntry_co, seq_type: SeqType = SeqType.NT
) -> Optional[bool]:
    consecutive_offsets = True

    for i in range(len(REGIONS) - 1):
        prev_region_end = _getattr(rearrangement, f"{REGIONS[i]}_end", seq_type)
        next_region_start = _getattr(rearrangement, f"{REGIONS[i+1]}_start", seq_type)
        if (
            prev_region_end is not None
            and prev_region_end > 0
            and next_region_start is not None
            and next_region_start > 0
        ):
            consecutive_offsets = consecutive_offsets and next_region_start == prev_region_end + 1

    start = 0
    end = 0
    for region in REGIONS:
        region_start = _getattr(rearrangement, f"{region}_start", seq_type)
        positive_region_start = region_start is not None and region_start > 0
        if positive_region_start:
            consecutive_offsets = consecutive_offsets and region_start > start and region_start > end

        region_end = _getattr(rearrangement, f"{region}_end", seq_type)
        positive_region_end = region_end is not None and region_end > 0
        if positive_region_end:
            consecutive_offsets = consecutive_offsets and region_end > start and region_end > end
            if positive_region_start:
                consecutive_offsets = consecutive_offsets and region_end >= region_start

        if positive_region_start:
            start = region_start
        if positive_region_end:
            end = region_end

    return consecutive_offsets


def validate_no_negative_offsets_inside_gene_alignment(
    rearrangement: AirrRearrangementEntry_co, gene: Gene, seq_type: SeqType = SeqType.NT
) -> Optional[bool]:
    regions_with_positive_start_end = [
        region
        for region in GENE_TO_REGIONS_START_END[gene]
        if _getattr(rearrangement, region, seq_type) is not None and _getattr(rearrangement, region, seq_type) > 0
    ]
    if "_".join(regions_with_positive_start_end) in "_".join(GENE_TO_REGIONS_START_END[gene]):
        return True
    return False


def validate_correct_region_offsets(
    rearrangement: AirrRearrangementEntry_co, seq_type: SeqType = SeqType.NT
) -> EntryValidationResults:
    def _set_validator_names(validation_results: EntryValidationResults) -> EntryValidationResults:
        return {f"correct_{region}{seq_type.value}_offsets": val for region, val in validation_results.items()}

    validation_results: EntryValidationResults = {region: None for region in REGIONS}

    if _getattr(rearrangement, "v_sequence_alignment", seq_type) is None:
        return _set_validator_names(validation_results)

    sequence = _getattr(rearrangement, "sequence", seq_type)
    fwr1_start = _getattr(rearrangement, "fwr1_start", seq_type)
    fwr1_end = _getattr(rearrangement, "fwr1_end", seq_type)
    if fwr1_start is not None and fwr1_end is not None and fwr1_end != -1:
        if fwr1_start == -1:
            fwr1_start = _getattr(rearrangement, "v_sequence_start", seq_type)
        validation_results["fwr1"] = sequence[fwr1_start - 1 : fwr1_end] == _getattr(rearrangement, "fwr1", seq_type)

    for region in REGIONS[1:-1]:
        # type: ignore
        region_start, region_end = _getattr(rearrangement, f"{region}_start", seq_type), _getattr(
            rearrangement, f"{region}_end", seq_type
        )
        if region_start is None or region_end is None or region_start == -1 or region_end == -1:
            continue
        validation_results[region] = sequence[region_start - 1 : region_end] == _getattr(
            rearrangement, region, seq_type
        )  # type: ignore

    if _getattr(rearrangement, "j_sequence_alignment", seq_type) is None:
        return _set_validator_names(validation_results)

    fwr4_start = _getattr(rearrangement, "fwr4_start", seq_type)
    fwr4_end = _getattr(rearrangement, "fwr4_end", seq_type)
    if fwr4_start != -1 and fwr4_start is not None and fwr4_end is not None:
        if fwr4_end == -1:
            fwr4_end = _getattr(rearrangement, "j_sequence_end", seq_type)
        validation_results["fwr4"] = sequence[fwr4_start - 1 : fwr4_end] == _getattr(rearrangement, "fwr4", seq_type)

    return _set_validator_names(validation_results)


def validate_correct_vj_in_frame(rearrangement: AirrRearrangementEntryNT) -> Optional[bool]:
    j_alignment_start, j_frame = rearrangement.j_alignment_start, rearrangement.j_frame
    v_alignment_start, v_frame = rearrangement.v_alignment_start, rearrangement.v_frame

    if (
        j_alignment_start is None
        or j_frame is None
        or v_alignment_start is None
        or v_frame is None
        or rearrangement.vj_in_frame is None
    ):
        return None

    vj_in_frame = ((j_alignment_start - 1 + j_frame) - (v_alignment_start - 1 + v_frame)) % 3 == 0

    return vj_in_frame == rearrangement.vj_in_frame


def validate_cdr3_in_junction(rearrangement: AirrRearrangementEntryNT) -> Optional[bool]:
    cdr3 = rearrangement.cdr3
    junction = rearrangement.junction
    cdr3_aa = rearrangement.cdr3_aa
    junction_aa = rearrangement.junction_aa
    if (not cdr3) or (not junction) or (not cdr3_aa) or (not junction_aa):
        return None
    return cdr3 in junction and cdr3_aa in junction_aa


def validate_cdr3_aa_in_junction_aa(rearrangement: AirrRearrangementEntryAA) -> Optional[bool]:
    cdr3_aa = rearrangement.cdr3_aa
    junction_aa = rearrangement.junction_aa
    if (not cdr3_aa) or (not junction_aa):
        return None
    return cdr3_aa in junction_aa


def validate_locus_as_in_v_gene(rearrangement: AirrRearrangementEntry_co) -> Optional[bool]:
    locus = rearrangement.locus
    v_call = rearrangement.v_call
    if not (locus) or (not v_call):
        return None
    return locus.upper() in v_call.upper()


def validate_gene_alignment(
    rearrangement: AirrRearrangementEntry_co, gene: Gene, seq_type: SeqType = SeqType.NT
) -> Optional[bool]:
    sequence = _getattr(rearrangement, "sequence", seq_type)
    sequence_start = _getattr(rearrangement, "sequence_start", seq_type, gene)
    sequence_end = _getattr(rearrangement, "sequence_end", seq_type, gene)
    sequence_alignment = _getattr(rearrangement, "sequence_alignment", seq_type, gene).replace("-", "")

    if (not sequence_start) or (not sequence_end) or (not sequence_alignment):
        return None
    return sequence[sequence_start - 1 : sequence_end] == sequence_alignment


def validate_no_empty_regions_in_gene_alignment(
    rearrangement: AirrRearrangementEntry_co, gene: Gene, seq_type: SeqType = SeqType.NT
) -> EntryValidationResults:
    def _set_validator_names(validation_results: EntryValidationResults) -> EntryValidationResults:
        return {
            f"no_empty_{region}{seq_type.value}_in_{gene.value.rstrip('_')}": val
            for region, val in validation_results.items()
        }

    gene_full_regions = [
        region
        for region, start_end_gene in REGION_TO_GENE.items()
        if start_end_gene["start"] == start_end_gene["end"] == gene
    ]
    validation_results: EntryValidationResults = {region: None for region in gene_full_regions}
    gene_sequence_alignment = _getattr(rearrangement, "sequence_alignment", seq_type, gene).replace("-", "")
    if gene_sequence_alignment is None:
        return _set_validator_names(validation_results)
    for region in gene_full_regions:
        if _getattr(rearrangement, region, seq_type) is not None:
            validation_results[region] = True
        else:
            validation_results[region] = False
    return _set_validator_names(validation_results)


def validate_no_empty_cdr3(rearrangement: AirrRearrangementEntry_co, seq_type: SeqType = SeqType.NT) -> Optional[bool]:
    if (
        _getattr(rearrangement, "v_sequence_alignment", seq_type) is None
        or _getattr(rearrangement, "j_sequence_alignment", seq_type) is None
    ):
        return None
    return _getattr(rearrangement, "cdr3", seq_type) is not None


def validate_cdr3_starts_inside_v_alignment(
    rearrangement: AirrRearrangementEntry_co, seq_type: SeqType = SeqType.NT
) -> Optional[bool]:
    cdr3_start = _getattr(rearrangement, "cdr3_start", seq_type)
    v_sequence_end = _getattr(rearrangement, "v_sequence_end", seq_type)
    if cdr3_start is None or v_sequence_end is None:
        return None
    if cdr3_start == -1 or cdr3_start > v_sequence_end:
        return False
    return True


def validate_conserved_residues_present(rearrangement: AirrRearrangementEntry_co) -> EntryValidationResults:
    imgt_conserved_residues = (
        IMGT_CONSERVED_RESIDUES_HEAVY if rearrangement.locus == "igh" else IMGT_CONSERVED_RESIDUES_LIGHT
    )
    numbering = rearrangement.scheme_residue_mapping
    validation_results: EntryValidationResults = {
        f"conserved_{aa}{pos}_present": None for pos, aa in imgt_conserved_residues.items()
    }
    if numbering is None:
        return validation_results

    return {
        f"conserved_{aa}{pos}_present": numbering.get(pos.split("_", maxsplit=1)[0]) == aa
        for pos, aa in imgt_conserved_residues.items()
    }


def validate_primary_sequence_in_sequence_alignment_aa(rearrangement: AirrRearrangementEntry_co) -> Optional[bool]:
    if rearrangement.sequence_alignment_aa is None or rearrangement.scheme_residue_mapping is None:
        return None
    return rearrangement.primary_sequence_aa in rearrangement.sequence_alignment_aa.replace("-", "")


def validate_no_insertion_next_to_deletion_aa(rearrangement: AirrRearrangementEntry_co) -> Optional[bool]:
    if rearrangement.sequence_aa_scheme_cigar is None:
        return None
    alignment_str = unfold_cigar(rearrangement.sequence_aa_scheme_cigar)
    return "ID" not in alignment_str and "DI" not in alignment_str


def validate_insertions_in_correct_places(rearrangement: AirrRearrangementEntry_co) -> Optional[bool]:
    if rearrangement.scheme_residue_mapping is None:
        return None

    scheme = Scheme(rearrangement.numbering_scheme)
    chain_type = ChainType.from_locus(Locus(rearrangement.locus))

    allowed_cdr_insertions: set[int] = set(get_indels_positions(scheme, chain_type))
    if scheme == Scheme.IMGT:
        allowed_cdr_insertions = {*allowed_cdr_insertions, *[ins + 1 for ins in allowed_cdr_insertions]}

    allowed_framework_insertions = get_framework_insertion_positions(scheme, chain_type)
    for scheme_pos in rearrangement.scheme_residue_mapping:
        split = scheme_pos.split(".")
        if len(split) > 1:
            if (int(split[0]) not in allowed_cdr_insertions) and (scheme_pos not in allowed_framework_insertions):
                return False
    return True


VALIDATORS: dict[str, Callable[[AirrRearrangementEntryNT], Optional[bool]]] = {
    "regions_in_aligned_sequence": validate_regions_in_aligned_sequence,
    "regions_aa_in_aligned_sequence_aa": partial(validate_regions_in_aligned_sequence, seq_type=SeqType.AA),
    "translated_regions_in_aligned_sequence_aa": validate_translated_regions_in_aligned_sequence_aa,
    "correct_vj_in_frame": validate_correct_vj_in_frame,
    "cdr3_in_junction": validate_cdr3_in_junction,
    "locus_as_in_v_gene": validate_locus_as_in_v_gene,
    "v_gene_alignment": partial(validate_gene_alignment, gene=Gene.V),
    "j_gene_alignment": partial(validate_gene_alignment, gene=Gene.J),
    "c_gene_alignment": partial(validate_gene_alignment, gene=Gene.C),
    "no_negative_offsets_inside_v_alignment": partial(validate_no_negative_offsets_inside_gene_alignment, gene=Gene.V),
    "no_negative_offsets_inside_j_alignment": partial(validate_no_negative_offsets_inside_gene_alignment, gene=Gene.J),
    "consecutive_offsets": validate_consecutive_offsets,
    "no_empty_cdr3": validate_no_empty_cdr3,
    "primary_sequence_in_sequence_alignment_aa": validate_primary_sequence_in_sequence_alignment_aa,
    "no_insertion_next_to_deletion_aa": validate_no_insertion_next_to_deletion_aa,
    "insertions_in_correct_places": validate_insertions_in_correct_places,
}
MULTI_VALIDATORS: list[Callable[[AirrRearrangementEntryNT], EntryValidationResults]] = [
    validate_correct_region_offsets,
    partial(validate_no_empty_regions_in_gene_alignment, gene=Gene.V),
    partial(validate_no_empty_regions_in_gene_alignment, gene=Gene.J),
]


def validate_airr_entry(rearrangement: AirrRearrangementEntryNT, scheme: Scheme) -> dict[str, Optional[bool]]:
    entry_validation_results = {}
    for validator_name, validator in VALIDATORS.items():
        try:
            entry_validation_results[validator_name] = validator(rearrangement)
        except Exception:  # pylint: disable=broad-except
            pass

    for multi_validator in MULTI_VALIDATORS:
        try:
            partial_validation_result = multi_validator(rearrangement)
            entry_validation_results.update(partial_validation_result)
        except Exception:  # pylint: disable=broad-except
            pass

    try:
        match scheme:
            case Scheme.IMGT:
                partial_validation_result = validate_conserved_residues_present(rearrangement)
                entry_validation_results.update(partial_validation_result)
            case _:
                pass
    except Exception:  # pylint: disable=broad-except
        pass

    return entry_validation_results


VALIDATORS_AA: dict[str, Callable[[AirrRearrangementEntryAA], Optional[bool]]] = {
    "regions_aa_in_aligned_sequence_aa": partial(validate_regions_in_aligned_sequence, seq_type=SeqType.AA),
    "cdr3_aa_in_junction_aa": validate_cdr3_aa_in_junction_aa,
    "locus_as_in_v_gene": validate_locus_as_in_v_gene,
    "v_gene_alignment_aa": partial(validate_gene_alignment, gene=Gene.V, seq_type=SeqType.AA),
    "j_gene_alignment_aa": partial(validate_gene_alignment, gene=Gene.J, seq_type=SeqType.AA),
    "c_gene_alignment_aa": partial(validate_gene_alignment, gene=Gene.C, seq_type=SeqType.AA),
    "no_negative_offsets_inside_v_alignment_aa": partial(
        validate_no_negative_offsets_inside_gene_alignment, gene=Gene.V, seq_type=SeqType.AA
    ),
    "no_negative_offsets_inside_j_alignment_aa": partial(
        validate_no_negative_offsets_inside_gene_alignment, gene=Gene.J, seq_type=SeqType.AA
    ),
    "consecutive_offsets_aa": partial(validate_consecutive_offsets, seq_type=SeqType.AA),
    "no_empty_cdr3_aa": partial(validate_no_empty_cdr3, seq_type=SeqType.AA),
    "primary_sequence_in_sequence_alignment_aa": validate_primary_sequence_in_sequence_alignment_aa,
    "no_insertion_next_to_deletion_aa": validate_no_insertion_next_to_deletion_aa,
    "insertions_in_correct_places": validate_insertions_in_correct_places,
}
MULTI_VALIDATORS_AA: list[Callable[[AirrRearrangementEntryAA], EntryValidationResults]] = [
    partial(validate_correct_region_offsets, seq_type=SeqType.AA),
    partial(validate_no_empty_regions_in_gene_alignment, gene=Gene.V, seq_type=SeqType.AA),
    partial(validate_no_empty_regions_in_gene_alignment, gene=Gene.J, seq_type=SeqType.AA),
]


def validate_airr_entry_aa(rearrangement: AirrRearrangementEntryAA, scheme: Scheme) -> dict[str, Optional[bool]]:
    entry_validation_results = {}
    for validator_name, validator in VALIDATORS_AA.items():
        try:
            entry_validation_results[validator_name] = validator(rearrangement)  # type: ignore
        except Exception:  # pylint: disable=broad-except
            pass

    for multi_validator in MULTI_VALIDATORS_AA:
        try:
            partial_validation_result = multi_validator(rearrangement)  # type: ignore
            entry_validation_results.update(partial_validation_result)
        except Exception:  # pylint: disable=broad-except
            pass

    try:
        match scheme:
            case Scheme.IMGT:
                partial_validation_result = validate_conserved_residues_present(rearrangement)
                entry_validation_results.update(partial_validation_result)
            case _:
                pass
    except Exception:  # pylint: disable=broad-except
        pass

    return entry_validation_results
