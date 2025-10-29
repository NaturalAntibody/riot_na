import os
from pathlib import Path
from typing import Literal, Optional

from riot_na.alignment.alignment_utils import unfold_cigar
from riot_na.config import GENE_DB_DIR
from riot_na.data.constants import AMINO_ACIDS
from riot_na.data.model import (
    AlignmentEntryAA,
    AlignmentsAA,
    AlignmentString,
    ChainType,
    Cigar,
    Locus,
    Organism,
    Scheme,
    SchemeAlignment,
)
from riot_na.data.scheme_definitions import get_legal_positions
from riot_na.data.scheme_mapping_facade import SchemeMappingFacade
from riot_na.schemes.collapse_alignment import collapse_alignment_str, collapse_ins_del
from riot_na.schemes.smooth_alignment import smooth_cdr_junctions

IMGT_REVERSE_ORDER_INSERTIONS = [33, 61, 112]


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


def get_scheme_residue_mapping(numbering: IntermediateNumbering) -> dict[str, str]:
    result: dict[str, str] = {}

    for scheme_id, residues in numbering.items():
        for residue, insertion in residues:
            scheme_id_str = str(scheme_id) if not insertion else f"{scheme_id}.{insertion}"
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
        allowed_species: Optional[list[Organism]] = None,
        db_dir: Path = GENE_DB_DIR,
    ):
        self.scheme_mapping_facades: dict[Scheme, SchemeMappingFacade] = {}

        if not allowed_species:
            allowed_species = [Organism.HOMO_SAPIENS, Organism.MUS_MUSCULUS, Organism.VICUGNA_PACOS]

        for scheme in Scheme:
            self.scheme_mapping_facades[scheme] = SchemeMappingFacade(scheme, allowed_species, db_dir)

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
