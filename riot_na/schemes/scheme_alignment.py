import os
from pathlib import Path
from typing import Literal, Optional, Sequence

from riot_na.alignment.alignment_utils import unfold_cigar
from riot_na.config import GENE_DB_DIR
from riot_na.data.constants import AMINO_ACIDS
from riot_na.data.model import (
    AlignmentEntryAA,
    AlignmentsAA,
    AlignmentString,
    ChainType,
    Cigar,
    CScheme,
    Locus,
    Organism,
    Scheme,
    SchemeAlignment,
)
from riot_na.data.scheme_definitions import (
    IMGT_C_DOMAINS,
    get_c_scheme_positions,
    get_legal_positions,
)
from riot_na.data.scheme_mapping_facade import CSchemeMappingFacade, SchemeMappingFacade
from riot_na.schemes.collapse_alignment import collapse_alignment_str, collapse_ins_del
from riot_na.schemes.smooth_alignment import smooth_cdr_junctions

IMGT_REVERSE_ORDER_INSERTIONS = [33, 61, 112]
# EU positions of germline residues enclosing the lower hinge (EU 231-237)
EU_LOWER_HINGE_ANCHORS = ("230", "238")


def validate_sequence(sequence):
    """
    Check whether a sequence is a protein sequence or if someone has submitted something nasty.
    """
    assert len(sequence) < 10000, "Sequence too long."
    unknown_aa = set(sequence.upper()) - AMINO_ACIDS
    assert not (unknown_aa), f"Unknown amino acid letter found in sequence: {list(unknown_aa)}"


def validate_alignment_len(aln: AlignmentEntryAA, gene: Literal["V", "J"], min_len: int):
    """
    Check whether alignment is not too short.
    """
    assert (
        aln.q_end - aln.q_start >= min_len
    ), f"{gene} sequence alignment AA too short, length: {aln.q_end - aln.q_start}"


def force_n_terminus_matches(alignment: AlignmentEntryAA) -> AlignmentEntryAA:
    # force matches on N terminus for V gene
    n_term_number_of_matches = min(alignment.q_start, alignment.t_start)
    fixed_cigar = Cigar(f"{n_term_number_of_matches}M" + alignment.cigar)

    q_start = alignment.q_start - n_term_number_of_matches
    t_start = alignment.t_start - n_term_number_of_matches

    return AlignmentEntryAA(
        target_id=alignment.target_id,
        alignment_score=None,
        seq_identity=None,
        e_value=None,
        q_start=q_start,
        q_end=alignment.q_end,
        t_start=t_start,
        t_end=alignment.t_end,
        cigar=fixed_cigar,
        species=alignment.species,
        locus=alignment.locus,
        q_seq=alignment.q_seq,
        t_seq=alignment.t_seq,
    )


def force_c_terminus_matches(
    query_sequence: str, target_sequence: str, alignment: AlignmentEntryAA
) -> AlignmentEntryAA:
    # force matches on C terminus for J gene
    # c_term_number_of_matches = min(len(query_sequence) - (alignment.q_end -alignment.q_start), len(target_sequence) - alignment.t_end)
    c_term_number_of_matches = min(len(query_sequence) - (alignment.q_end), len(target_sequence) - alignment.t_end)
    fixed_cigar = Cigar(alignment.cigar + f"{c_term_number_of_matches}M")

    q_end = alignment.q_end + c_term_number_of_matches
    t_end = alignment.t_end + c_term_number_of_matches

    return AlignmentEntryAA(
        target_id=alignment.target_id,
        alignment_score=None,
        seq_identity=None,
        e_value=None,
        q_start=alignment.q_start,
        q_end=q_end,
        t_start=alignment.t_start,
        t_end=t_end,
        cigar=fixed_cigar,
        species=alignment.species,
        locus=alignment.locus,
        q_seq=alignment.q_seq,
        t_seq=alignment.t_seq,
    )


def force_n_terminus_del_ins(aln: AlignmentEntryAA, query_gene_alignment_str: AlignmentString) -> AlignmentString:
    # add deletions and insertions on N terminus of alignment so the cigar covers whole gene
    # this is needed to infer the last scheme position of V gene and first position of J gene
    if aln.t_start > 0:
        query_gene_alignment_str = AlignmentString(aln.t_start * "N" + query_gene_alignment_str)

    if aln.q_start > 0:
        query_gene_alignment_str = AlignmentString(aln.q_start * "S" + query_gene_alignment_str)
    return query_gene_alignment_str


def _merge_cigars(
    query_gene_alignment_str: AlignmentString,
    gene_scheme_alignment_str: AlignmentString,
    query_start: int = 0,
    target_start: int = 0,
) -> AlignmentString:
    # I M -> I: it1++ -ok
    # I D -> M:it1++, it2++ -ok
    # I I -> I: it2++ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # M M ->M:it1++, it2++ -ok
    # M D -> D: it2++ -ok
    # M I -> I: it1++, it2++ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # D D -> D: it1++ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # D M -> : D : it1++, it2++
    # D I -> X: it1++, it2++

    #     'SSSSSSSSSSSSSNNNMMMMMMMM'
    # 'NNNNNNNNNNNNNNNMMMDMMMMMMMMM'

    query_gene_pos = query_start + target_start
    gene_scheme_pos = target_start
    res_cigar = ""

    while query_gene_pos < len(query_gene_alignment_str) and gene_scheme_pos < len(gene_scheme_alignment_str):
        query_gene_op = query_gene_alignment_str[query_gene_pos]
        gene_scheme_op = gene_scheme_alignment_str[gene_scheme_pos]

        match (query_gene_op, gene_scheme_op):
            case ("I", "M"):
                query_gene_pos = query_gene_pos + 1
                res_cigar = res_cigar + "I"
            case ("I", "D"):
                query_gene_pos = query_gene_pos + 1
                gene_scheme_pos = gene_scheme_pos + 1
                res_cigar = res_cigar + "M"
            case ("I", "I"):
                query_gene_pos = query_gene_pos + 1
                res_cigar = res_cigar + "I"
            case ("M", "M"):
                query_gene_pos = query_gene_pos + 1
                gene_scheme_pos = gene_scheme_pos + 1
                res_cigar = res_cigar + "M"
            case ("M", "D"):
                gene_scheme_pos = gene_scheme_pos + 1
                res_cigar = res_cigar + "D"
            case ("M", "I"):
                query_gene_pos = query_gene_pos + 1
                gene_scheme_pos = gene_scheme_pos + 1
                res_cigar = res_cigar + "I"
            case ("D", "D"):
                gene_scheme_pos = gene_scheme_pos + 1
                res_cigar = res_cigar + "D"
            case ("D", "M"):
                query_gene_pos = query_gene_pos + 1
                gene_scheme_pos = gene_scheme_pos + 1
                res_cigar = res_cigar + "D"
            case ("D", "I"):
                query_gene_pos = query_gene_pos + 1
                gene_scheme_pos = gene_scheme_pos + 1
            case ("N", "M"):
                query_gene_pos = query_gene_pos + 1
                gene_scheme_pos = gene_scheme_pos + 1
                res_cigar = res_cigar + "N"
            case ("N", "D"):
                gene_scheme_pos = gene_scheme_pos + 1
                res_cigar = res_cigar + "N"
            case ("N", "I"):
                query_gene_pos = query_gene_pos + 1
                gene_scheme_pos = gene_scheme_pos + 1
            case ("S", "M"):
                query_gene_pos = query_gene_pos + 1
                res_cigar = res_cigar + "S"
            case ("S", "I"):
                query_gene_pos = query_gene_pos + 1
                res_cigar = res_cigar + "S"
            case ("S", "D"):
                query_gene_pos = query_gene_pos + 1
                gene_scheme_pos = gene_scheme_pos + 1

    while query_gene_pos < len(query_gene_alignment_str):
        query_gene_op = query_gene_alignment_str[query_gene_pos]
        assert query_gene_op not in {"M", "D"}, "Query gene alignment is longer than gene scheme alignment"
        res_cigar += query_gene_op
        query_gene_pos = query_gene_pos + 1

    return AlignmentString(res_cigar)


def infer_last_scheme_position_aligned(query_scheme_aln_str: AlignmentString) -> int:
    # query =   "CAAAAAAAAAAAA"
    # gene =    "BAAAAAAABC"
    # aln =    "SNMMMMMMM"
    # g_sch =   "MMMDMMMMMMM"
    # q_sch =  "SNMMDMMMM"

    # scheme space:
    # query =  "C-AA AAAAAAAAAA"
    # gene =   "-BAA AAAAABC"
    # aln =    "SNMM MMMMM"
    # g_sch =   "MMMDMMMMMMM"
    # q_sch =  "IDMMDMMMMM"
    # last_pos= "        |"  # 9 (1-based)

    matches = query_scheme_aln_str.count("M")
    deletions = query_scheme_aln_str.count("D") + query_scheme_aln_str.count("N")

    last_position_sch = matches + deletions  # 1-based
    return last_position_sch


def merge_v_j(
    v_query_scheme_cigar: AlignmentString,
    j_query_scheme_cigar: AlignmentString,
    j_gene_start_on_scheme: int,
) -> AlignmentString:
    last_v_match_pos = v_query_scheme_cigar.rfind("M")
    first_j_match_pos = j_query_scheme_cigar.find("M")

    last_aligned_v_position_on_scheme = infer_last_scheme_position_aligned(
        AlignmentString(v_query_scheme_cigar[: last_v_match_pos + 1])
    )

    # remaining part of alignmetns that does not match to the scheme
    middle_v_scheme_cigar = v_query_scheme_cigar[last_v_match_pos + 1 :]
    middle_j_scheme_cigar = j_query_scheme_cigar[:first_j_match_pos]

    scheme_deletions = (
        j_gene_start_on_scheme
        + middle_j_scheme_cigar.count("N")
        + middle_j_scheme_cigar.count("D")
        - last_aligned_v_position_on_scheme
        - 1
    )
    query_insertions = (
        middle_v_scheme_cigar.count("I") + middle_j_scheme_cigar.count("S") + middle_j_scheme_cigar.count("I")
    )
    middle_cigar = collapse_ins_del(scheme_deletions * ["D"] + query_insertions * ["I"])

    final_cigar = AlignmentString(
        f"{v_query_scheme_cigar[:last_v_match_pos + 1]}{middle_cigar}{j_query_scheme_cigar[first_j_match_pos:]}"
    )
    return final_cigar


IntermediateNumbering = dict[int, list[tuple[str, int]]]


def produce_numbering(
    query: str, alignment_str: AlignmentString, query_start: Optional[int] = None
) -> IntermediateNumbering:
    scheme_position = 0
    query_position = query_start or 0
    insertion_counter = 0

    numbering: IntermediateNumbering = {}
    for op in alignment_str:
        if op == "M":
            insertion_counter = 0
            query_position = query_position + 1
            scheme_position = scheme_position + 1
            position_residues = numbering.get(scheme_position, [])
            position_residues.append((query[query_position - 1], insertion_counter))
            numbering[scheme_position] = position_residues
        elif op in {"D", "N"}:
            insertion_counter = 0
            scheme_position = scheme_position + 1
        elif op == "I":
            insertion_counter = insertion_counter + 1
            query_position = query_position + 1
            position_residues = numbering.get(scheme_position, [])
            position_residues.append((query[query_position - 1], insertion_counter))
            numbering[scheme_position] = position_residues

    return numbering


def scheme_alignment(
    v_aligment: AlignmentEntryAA,
    v_gene_scheme_alignment_str: AlignmentString,
    j_aligment: Optional[AlignmentEntryAA],
    j_gene_scheme_alignment_str: Optional[AlignmentString],
    j_gene_start_on_scheme: Optional[int],
    scheme: Scheme,
    extend_alignment: bool = True,
) -> SchemeAlignment:
    fixed_v_aln = v_aligment
    if extend_alignment:
        fixed_v_aln = force_n_terminus_matches(v_aligment)

    v_query_gene_alignment_str = unfold_cigar(fixed_v_aln.cigar)
    v_query_gene_alignment_str = collapse_alignment_str(v_query_gene_alignment_str)
    v_query_gene_alignment_str = force_n_terminus_del_ins(fixed_v_aln, v_query_gene_alignment_str)

    v_query_scheme_alignment_str = _merge_cigars(v_query_gene_alignment_str, v_gene_scheme_alignment_str)
    query_scheme_alignment_str = v_query_scheme_alignment_str

    if j_aligment:
        assert j_gene_scheme_alignment_str is not None
        assert j_gene_start_on_scheme is not None
        assert j_aligment.t_seq is not None

        if j_aligment.q_end - j_aligment.q_start >= int(os.environ.get("ALIGNMENT_LENGTH_THRESHOLD_J_AA", 5)):
            masked = v_aligment.q_seq[v_aligment.q_end :]

            # masked = v_aligment.q_seq[v_aligment.q_end - v_aligment.q_start :]
            if extend_alignment:
                j_aligment = force_c_terminus_matches(masked, j_aligment.t_seq, j_aligment)

            j_query_gene_alignment_str = unfold_cigar(j_aligment.cigar)
            j_query_gene_alignment_str = collapse_alignment_str(j_query_gene_alignment_str)
            j_query_gene_alignment_str = force_n_terminus_del_ins(j_aligment, j_query_gene_alignment_str)

            j_query_scheme_alignment_str = _merge_cigars(j_query_gene_alignment_str, j_gene_scheme_alignment_str)

            query_scheme_alignment_str = merge_v_j(
                v_query_scheme_cigar=query_scheme_alignment_str,
                j_query_scheme_cigar=j_query_scheme_alignment_str,
                j_gene_start_on_scheme=j_gene_start_on_scheme,
            )

    chain_type = ChainType.HEAVY if v_aligment.locus == Locus.IGH else ChainType.LIGHT
    collapsed_query_scheme_alignment_str = collapse_alignment_str(query_scheme_alignment_str, ordered=True)
    smoothed_collapsed_query_scheme_alignment_str = smooth_cdr_junctions(
        collapsed_query_scheme_alignment_str, chain_type, scheme
    )

    return SchemeAlignment(
        target_id=scheme.value,
        q_start=fixed_v_aln.q_start,
        q_end=j_aligment.q_end + fixed_v_aln.q_end if j_aligment else fixed_v_aln.q_end,
        alignment_str=smoothed_collapsed_query_scheme_alignment_str,
    )


def clip_to_gene_scheme_alignment(
    query_gene_alignment_str: AlignmentString, gene_scheme_alignment_str: AlignmentString
) -> AlignmentString:
    # germline residues past the end of the gene scheme alignment (e.g. membrane exons of C genes) have no
    # scheme positions, drop the part of the query aligned to them
    numbered_gene_length = gene_scheme_alignment_str.count("M") + gene_scheme_alignment_str.count("I")
    gene_position = 0

    for i, op in enumerate(query_gene_alignment_str):
        if gene_position == numbered_gene_length:
            return AlignmentString(query_gene_alignment_str[:i])
        if op in {"M", "D", "N"}:
            gene_position = gene_position + 1

    return query_gene_alignment_str


def get_germline_scheme_indices(
    gene_scheme_alignment_str: AlignmentString, scheme_positions: Sequence[str]
) -> dict[str, int]:
    """Indices of germline residues on scheme positions."""
    indices: dict[str, int] = {}
    scheme_index = germline_index = 0

    for op in gene_scheme_alignment_str:
        if op == "M":
            indices[scheme_positions[scheme_index]] = germline_index
        if op in {"M", "D"}:
            scheme_index = scheme_index + 1
        if op in {"M", "I"}:
            germline_index = germline_index + 1

    return indices


def get_lower_hinge_anchors(
    eu_gene_scheme_alignment_str: Optional[AlignmentString], chain_type: ChainType
) -> Optional[tuple[int, int]]:
    """Indices of germline residues on EU 230 and 238, None if the germline has no lower hinge."""
    if not eu_gene_scheme_alignment_str:
        return None

    indices = get_germline_scheme_indices(eu_gene_scheme_alignment_str, get_c_scheme_positions(CScheme.EU, chain_type))
    first_anchor, last_anchor = EU_LOWER_HINGE_ANCHORS
    if first_anchor not in indices or last_anchor not in indices:
        return None
    return indices[first_anchor], indices[last_anchor]


def _count_identities(query_gene_alignment_str: str, query: str, germline: str) -> int:
    identities = 0
    query_position = germline_position = 0

    for op in query_gene_alignment_str:
        if op == "M" and query[query_position] == germline[germline_position]:
            identities = identities + 1
        if op in {"M", "I"}:
            query_position = query_position + 1
        if op in {"M", "D"}:
            germline_position = germline_position + 1

    return identities


def place_lower_hinge_indels(
    query_gene_alignment_str: AlignmentString,
    query: str,
    germline: str,
    query_start: int,
    anchors: tuple[int, int],
) -> AlignmentString:
    """
    Smith-Waterman places gaps where substitution scores are the highest, so a residue of the lower hinge
    (EU 231-237) missing in the query or an extra one is placed depending on the mutations next to it, e.g. IGHG1
    with E233P/L234V/L235A/ΔG236 is aligned with L234 missing. It is placed right before G237 instead, where IGHG2
    lacks G236 in EU and Kabat mappings (the Fc engineering convention, IMGT places the IGHG2 deletion at CH2 1.4),
    unless the alignment pairs more identical residues. Lower hinges differing from the germline by more residues
    are left as aligned.

    anchors - indices of germline residues on EU 230 and 238, both have to be aligned.
    """
    first_anchor, last_anchor = anchors
    # germline index of an anchor -> its index in the alignment string and the index of the query residue paired with it
    aligned_anchors: dict[int, tuple[int, int]] = {}
    query_position, germline_position = query_start, 0

    for i, op in enumerate(query_gene_alignment_str):
        if op == "M" and germline_position in anchors:
            aligned_anchors[germline_position] = (i, query_position)
        if op in {"M", "I"}:
            query_position = query_position + 1
        if op in {"M", "D", "N"}:
            germline_position = germline_position + 1

    if len(aligned_anchors) < 2:
        return query_gene_alignment_str

    (first_op_index, first_query_position), (last_op_index, last_query_position) = (
        aligned_anchors[first_anchor],
        aligned_anchors[last_anchor],
    )
    hinge_query = query[first_query_position + 1 : last_query_position]
    hinge_germline = germline[first_anchor + 1 : last_anchor]
    paired = min(len(hinge_query), len(hinge_germline))
    if abs(len(hinge_query) - len(hinge_germline)) != 1 or not paired:
        return query_gene_alignment_str

    # residues are paired from EU 231 on, the last one with G237
    indel = "I" if len(hinge_query) > len(hinge_germline) else "D"
    placed = (paired - 1) * "M" + indel + "M"

    aligned = query_gene_alignment_str[first_op_index + 1 : last_op_index]
    if _count_identities(placed, hinge_query, hinge_germline) < _count_identities(aligned, hinge_query, hinge_germline):
        return query_gene_alignment_str

    return AlignmentString(
        query_gene_alignment_str[: first_op_index + 1] + placed + query_gene_alignment_str[last_op_index:]
    )


def c_scheme_alignment(
    c_alignment: AlignmentEntryAA,
    c_gene_scheme_alignment_str: AlignmentString,
    c_scheme: CScheme,
    query_offset: int = 0,
    lower_hinge_anchors: Optional[tuple[int, int]] = None,
) -> SchemeAlignment:
    """
    Align query C region to the C scheme. query_offset is the position of the C alignment query
    (the query left after J alignment) in the whole query. lower_hinge_anchors are indices of germline residues on
    EU 230 and 238 of IGHG germlines.
    """
    c_query_gene_alignment_str = collapse_alignment_str(unfold_cigar(c_alignment.cigar))
    # unaligned query residues before the C alignment are not numbered, unlike V no S ops are added
    c_query_gene_alignment_str = AlignmentString(c_alignment.t_start * "N" + c_query_gene_alignment_str)
    if lower_hinge_anchors:
        c_query_gene_alignment_str = place_lower_hinge_indels(
            c_query_gene_alignment_str, c_alignment.q_seq, c_alignment.t_seq, c_alignment.q_start, lower_hinge_anchors
        )
    c_query_gene_alignment_str = clip_to_gene_scheme_alignment(c_query_gene_alignment_str, c_gene_scheme_alignment_str)

    # germline insertions next to germline deletions (e.g. Kabat 107A of lambda C after missing 107) are kept,
    # so unlike V the alignment is not collapsed
    c_query_scheme_alignment_str = _merge_cigars(c_query_gene_alignment_str, c_gene_scheme_alignment_str)

    q_start = query_offset + c_alignment.q_start
    q_end = q_start + c_query_gene_alignment_str.count("M") + c_query_gene_alignment_str.count("I")

    return SchemeAlignment(
        target_id=c_scheme.value,
        q_start=q_start,
        q_end=q_end,
        alignment_str=c_query_scheme_alignment_str,
    )


def _move_position_to_last_insertion(numbering: IntermediateNumbering, scheme_id: int) -> IntermediateNumbering:
    # residues following a position whose insertions precede it can't be its insertions, the position goes to the
    # last of them and the residues before it become insertions of the previous position
    residues = numbering.get(scheme_id, [])
    if not any(insertion for _, insertion in residues):
        return numbering

    previous_residues = numbering.get(scheme_id - 1, [])
    previous_insertions = len([residue for residue in previous_residues if residue[1]])
    numbering[scheme_id - 1] = previous_residues + [
        (residue, previous_insertions + i + 1) for i, (residue, _) in enumerate(residues[:-1])
    ]
    numbering[scheme_id] = [(residues[-1][0], 0)]
    return numbering


def _move_insertions_to_next_position(
    numbering: IntermediateNumbering, scheme_id: int, keep: int = 0
) -> IntermediateNumbering:
    # insertions preceding the next position are numbered towards it, the last one is .1
    residues = numbering.get(scheme_id, [])
    moved = [residue for residue in residues if residue[1]][keep:]
    if not moved:
        return numbering

    remaining = residues[: len(residues) - len(moved)]
    if remaining:
        numbering[scheme_id] = remaining
    else:
        del numbering[scheme_id]

    numbering[scheme_id + 1] = [(residue, len(moved) - i) for i, (residue, _) in enumerate(moved)] + numbering.get(
        scheme_id + 1, []
    )
    return numbering


def fix_imgt_c_numbering(numbering: IntermediateNumbering, scheme_positions: Sequence[str]) -> IntermediateNumbering:
    """
    produce_numbering attaches insertions to the preceding scheme position. In IMGT C-DOMAIN numbering insertions
    before position 1 (1.x) belong to it, and insertions between 84 and 85 are split: 84.x follow 84, 85.x
    precede 85, the odd one goes to 84 (IGKC 84.1-84.5, 85.4-85.1). Residues inserted after position 1 or 85 can't
    be numbered as their insertions, so the position goes to the last inserted residue and the residues before it
    become insertions before it (1.x, or 84.x and 85.x split as above).
    """
    scheme_ids = {position: i + 1 for i, position in enumerate(scheme_positions)}

    for domain in IMGT_C_DOMAINS:
        if f"{domain}_1" not in scheme_ids:
            continue

        position_1 = scheme_ids[f"{domain}_1"]
        numbering = _move_position_to_last_insertion(numbering, position_1)
        numbering = _move_insertions_to_next_position(numbering, position_1 - 1)

        numbering = _move_position_to_last_insertion(numbering, scheme_ids[f"{domain}_85"])
        position_84 = scheme_ids[f"{domain}_84"]
        de_turn_insertions = len([residue for residue in numbering.get(position_84, []) if residue[1]])
        numbering = _move_insertions_to_next_position(
            numbering, position_84, keep=de_turn_insertions - de_turn_insertions // 2
        )

    return dict(sorted(numbering.items()))


def produce_c_numbering(
    query: str, c_sch_alignment: SchemeAlignment, c_scheme: CScheme, chain_type: ChainType
) -> dict[str, str]:
    """
    Scheme residue mapping of query C region aligned to the C scheme.
    """
    c_scheme_positions = get_c_scheme_positions(c_scheme, chain_type)

    c_numbering = produce_numbering(query, c_sch_alignment.alignment_str, c_sch_alignment.q_start)
    if c_scheme == CScheme.IMGT:
        c_numbering = fix_imgt_c_numbering(c_numbering, c_scheme_positions)

    return get_scheme_residue_mapping(c_numbering, c_scheme_positions)


def get_scheme_residue_mapping(
    numbering: IntermediateNumbering, scheme_positions: Optional[Sequence[str]] = None
) -> dict[str, str]:
    """
    scheme_positions - labels of scheme positions (1-based numbering indices), by default positions are labeled
    with their indices
    """
    result: dict[str, str] = {}

    for scheme_id, residues in numbering.items():
        scheme_position = scheme_positions[scheme_id - 1] if scheme_positions else str(scheme_id)
        for residue, insertion in residues:
            scheme_id_str = scheme_position if not insertion else f"{scheme_position}.{insertion}"
            result[scheme_id_str] = residue

    return result


def sort_imgt_numbering(numbering: IntermediateNumbering) -> IntermediateNumbering:
    for pos in IMGT_REVERSE_ORDER_INSERTIONS:
        if pos in numbering:
            numbering[pos] = sorted(numbering[pos], key=lambda x: x[1], reverse=True)
    return numbering


def get_positional_scheme_mapping(scheme_residue_mapping: dict[str, str], query_start: int = 0) -> dict[int, str]:
    return dict(enumerate(scheme_residue_mapping.keys(), query_start))


def _fix_imgt_cdr_numbering(
    numbering: IntermediateNumbering,
    insertion_position: int,
    swap_priority: bool = False,
) -> IntermediateNumbering:
    insertions = numbering.get(insertion_position, [])
    if len(insertions) <= 1:
        return numbering

    new_insertions = numbering.get(insertion_position + 1, [])
    assert len(new_insertions) == 1

    new_insertion_counter = 1
    while True:
        len_diff = len(new_insertions) - len(insertions) if swap_priority else len(insertions) - len(new_insertions)
        if len_diff in [0, 1]:
            break
        residue, _ = insertions.pop()
        new_insertions.append((residue, new_insertion_counter))
        new_insertion_counter += 1

    numbering[insertion_position] = insertions
    numbering[insertion_position + 1] = new_insertions
    return numbering


def fix_imgt_cdrs_numbering(numbering: IntermediateNumbering) -> IntermediateNumbering:
    # CDR1
    numbering = _fix_imgt_cdr_numbering(numbering, insertion_position=32)
    # CDR2
    numbering = _fix_imgt_cdr_numbering(numbering, insertion_position=60)
    # CDR3
    numbering = _fix_imgt_cdr_numbering(numbering, insertion_position=111, swap_priority=True)
    return numbering


def _get_j_gene_start_on_scheme(j_gene_scheme_alignment: AlignmentString, chain_type: ChainType, scheme: Scheme) -> int:
    legal_positions = get_legal_positions(chain_type, scheme)
    j_start_on_scheme = legal_positions - j_gene_scheme_alignment.count("M") - j_gene_scheme_alignment.count("D") + 1
    return j_start_on_scheme


class SchemeAligner:
    def __init__(
        self,
        allowed_species: Optional[tuple[Organism, ...]] = None,
        db_dir: Path = GENE_DB_DIR,
    ):
        self.scheme_mapping_facades: dict[Scheme, SchemeMappingFacade] = {}

        if not allowed_species:
            allowed_species = (Organism.HOMO_SAPIENS, Organism.MUS_MUSCULUS, Organism.VICUGNA_PACOS)

        for scheme in Scheme:
            self.scheme_mapping_facades[scheme] = SchemeMappingFacade(scheme, allowed_species, db_dir)

        self.c_scheme_mapping_facades: dict[CScheme, CSchemeMappingFacade] = {
            c_scheme: CSchemeMappingFacade(c_scheme, db_dir) for c_scheme in CScheme
        }

    def align_to_scheme(
        self, aa_alignments: AlignmentsAA, scheme: Scheme, extend_alignment: bool = True
    ) -> Optional[SchemeAlignment]:
        if not aa_alignments.v:
            return None

        organism = aa_alignments.v.species
        v_gene_scheme_alignment_str = self.scheme_mapping_facades[scheme].get_mapping(
            organism, aa_alignments.v.target_id
        )
        j_gene_scheme_alignment_str = (
            self.scheme_mapping_facades[scheme].get_mapping(organism, aa_alignments.j.target_id)
            if aa_alignments.j
            else None
        )

        locus = aa_alignments.v.locus
        chain_type = ChainType.from_locus(locus)

        j_start_on_scheme = (
            _get_j_gene_start_on_scheme(j_gene_scheme_alignment_str, chain_type, scheme)
            if j_gene_scheme_alignment_str
            else None
        )

        sch_alignment = scheme_alignment(
            aa_alignments.v,
            v_gene_scheme_alignment_str,
            aa_alignments.j,
            j_gene_scheme_alignment_str,
            j_start_on_scheme,
            scheme=scheme,
            extend_alignment=extend_alignment,
        )

        return sch_alignment

    def align_c_to_scheme(self, aa_alignments: AlignmentsAA, c_scheme: CScheme) -> Optional[SchemeAlignment]:
        """
        Align C region to the C scheme. Returns None if there is no C alignment or the C gene has no scheme mapping.
        """
        if not (aa_alignments.v and aa_alignments.c):
            return None

        c_gene_scheme_alignment_str = self.c_scheme_mapping_facades[c_scheme].get_mapping(
            aa_alignments.c.species, aa_alignments.c.target_id
        )
        if c_gene_scheme_alignment_str is None:
            return None

        # C is aligned to the query left after V and J alignments (or V alignment, if J was not found)
        c_query_offset = aa_alignments.v.q_end + (aa_alignments.j.q_end if aa_alignments.j else 0)

        # the lower hinge is found by EU positions in all C schemes
        lower_hinge_anchors = get_lower_hinge_anchors(
            self.c_scheme_mapping_facades[CScheme.EU].get_mapping(aa_alignments.c.species, aa_alignments.c.target_id),
            ChainType.from_locus(aa_alignments.c.locus),
        )

        return c_scheme_alignment(
            aa_alignments.c, c_gene_scheme_alignment_str, c_scheme, c_query_offset, lower_hinge_anchors
        )


if __name__ == "__main__":
    import json
    from dataclasses import asdict

    sch_aligner = SchemeAligner(allowed_species=None)

    v_almt = AlignmentEntryAA(
        target_id="IGHV4-59*01",
        alignment_score=None,
        seq_identity=None,
        e_value=None,
        q_start=6,
        q_end=82,
        q_seq="LRRPCPTCSISDGSISSYYWNWIRQSPGKGLEWIGHIHYSGSTHYNPSLQSRVSISIDTSKNHFSLKLRSVTAVDTAVYYCARWGHFDTSGYFVVDYWGQGTLVTVSS",
        t_start=20,
        t_end=96,
        t_seq="QVQLQESGPGLVKPSETLSLTCTVSGGSISSYYWSWIRQPPGKGLEWIGYIYYSGSTNYNPSLKSRVTISVDTSKNQFSLKLSSVTAADTAVYYCAR",
        cigar=Cigar("76M"),
        locus=Locus.IGH,
        species=Organism.HOMO_SAPIENS,
    )

    j_almt = AlignmentEntryAA(
        target_id="IGHJ4*02",
        alignment_score=None,
        seq_identity=None,
        e_value=None,
        q_start=16,
        q_end=23,
        q_seq="RWGHFDTSGYFVVDYWGQGTLVTVSS",
        t_start=5,
        t_end=12,
        t_seq="YFDYWGQGTLVTVSS",
        cigar=Cigar("7M"),
        species=Organism.HOMO_SAPIENS,
        locus=Locus.IGH,
    )

    aa_almts = AlignmentsAA(v=v_almt, j=j_almt)
    scheme_aln = sch_aligner.align_to_scheme(aa_almts, Scheme.MARTIN)

    if scheme_aln:
        print(json.dumps(asdict(scheme_aln)))
