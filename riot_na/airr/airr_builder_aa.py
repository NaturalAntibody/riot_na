from typing import Optional

from riot_na.airr.airr_validator import validate_airr_entry_aa
from riot_na.alignment.alignment_utils import (
    align_sequences,
    fold_cigar,
    offset_alignments,
    unfold_cigar,
)
from riot_na.data.model import (
    AirrRearrangementEntryAA,
    AlignmentEntryAA,
    AlignmentString,
    RegionOffsetsAA,
    Scheme,
    SchemeAlignment,
    SegmentedAirrRearrangementEntryAA,
)


class AirrBuilderAA:  # pylint: disable=too-many-instance-attributes
    def __init__(self, sequence_header: str, sequence: str, scheme: Scheme):
        self.sequence_header = sequence_header
        self.sequence_aa = sequence
        self.scheme = scheme

        self.rearrangement = AirrRearrangementEntryAA(
            sequence_header=sequence_header, sequence_aa=sequence, numbering_scheme=scheme.value
        )

        self.v_gene_sequence_aa: Optional[str] = None
        self.v_gene_alignment_aa: Optional[AlignmentEntryAA] = None

        self.j_gene_alignment_aa: Optional[AlignmentEntryAA] = None
        self.j_gene_sequence_aa: Optional[str] = None

        self.c_gene_alignment_aa: Optional[AlignmentEntryAA] = None
        self.c_gene_sequence_aa: Optional[str] = None

        self.aa_offsets: Optional[RegionOffsetsAA] = None

    def with_v_gene_alignment_aa(self, v_aln: AlignmentEntryAA):
        self.v_gene_sequence_aa = v_aln.t_seq
        self.v_gene_alignment_aa = v_aln

        if v_aln is None:
            return self

        self.rearrangement.v_call = v_aln.target_id
        self.rearrangement.locus = v_aln.locus.value
        self.rearrangement.locus_species = v_aln.species.value
        # alignment positions are 0-indexed, add 1 to start index to convert to 1-based
        self.rearrangement.v_sequence_start_aa = v_aln.q_start + 1
        self.rearrangement.v_sequence_end_aa = v_aln.q_end

        self.rearrangement.v_germline_start_aa = v_aln.t_start + 1
        self.rearrangement.v_germline_end_aa = v_aln.t_end

        v_sequence_segment = self.sequence_aa[v_aln.q_start : v_aln.q_end]
        v_germline_segment = v_aln.t_seq[v_aln.t_start : v_aln.t_end]

        v_sequence_alignment, v_germline_alignment = align_sequences(
            v_sequence_segment, v_germline_segment, unfold_cigar(v_aln.cigar)
        )

        self.rearrangement.v_sequence_alignment_aa = v_sequence_alignment
        self.rearrangement.v_germline_alignment_aa = v_germline_alignment

        # 1 based, hence need to add 1
        self.rearrangement.v_alignment_start_aa = 1
        self.rearrangement.v_alignment_end_aa = len(v_sequence_alignment)

        self.rearrangement.v_score_aa = v_aln.alignment_score
        self.rearrangement.v_cigar_aa = v_aln.cigar
        self.rearrangement.v_support_aa = v_aln.e_value
        self.rearrangement.v_identity_aa = v_aln.seq_identity

        return self

    def with_j_gene_alignment_aa(self, j_aln: AlignmentEntryAA):
        self.j_gene_sequence_aa = j_aln.t_seq

        if j_aln is not None:
            assert self.v_gene_alignment_aa is not None

            j_aln = offset_alignments(self.v_gene_alignment_aa.q_end, j_aln)

            self.j_gene_alignment_aa = j_aln

            self.rearrangement.j_call = j_aln.target_id
            # alignment positions are 0-indexed, add 1 to start index to convert to 1-based
            self.rearrangement.j_sequence_start_aa = j_aln.q_start + 1
            self.rearrangement.j_sequence_end_aa = j_aln.q_end

            self.rearrangement.j_germline_start_aa = j_aln.t_start + 1
            self.rearrangement.j_germline_end_aa = j_aln.t_end

            j_sequence_segment = self.sequence_aa[j_aln.q_start : j_aln.q_end]
            j_germline_segment = j_aln.t_seq[j_aln.t_start : j_aln.t_end]

            j_sequence_alignment, j_germline_alignment = align_sequences(
                j_sequence_segment, j_germline_segment, unfold_cigar(j_aln.cigar)
            )

            self.rearrangement.j_sequence_alignment_aa = j_sequence_alignment
            self.rearrangement.j_germline_alignment_aa = j_germline_alignment

            # 1 based, hence need to add 1
            assert self.rearrangement.v_sequence_start_aa is not None
            v_alignment_str = unfold_cigar(self.v_gene_alignment_aa.cigar)
            deletions_on_v = v_alignment_str.count("D")
            self.rearrangement.j_alignment_start_aa = j_aln.q_start + 1 + deletions_on_v
            self.rearrangement.j_alignment_end_aa = j_aln.q_end + deletions_on_v

            self.rearrangement.j_score_aa = j_aln.alignment_score
            self.rearrangement.j_cigar_aa = j_aln.cigar
            self.rearrangement.j_support_aa = j_aln.e_value
            self.rearrangement.j_identity_aa = j_aln.seq_identity

        return self

    def with_c_gene_alignment_aa(self, c_aln: AlignmentEntryAA):

        self.c_gene_sequence_aa = c_aln.t_seq

        if c_aln is not None:
            assert self.v_gene_alignment_aa is not None
            assert self.j_gene_alignment_aa is not None
            assert self.j_gene_alignment_aa.q_end is not None
            assert self.j_gene_alignment_aa.q_start is not None
            assert self.j_gene_alignment_aa.t_start is not None
            assert self.j_gene_alignment_aa.t_end is not None
            assert self.j_gene_alignment_aa.target_id is not None
            assert self.j_gene_alignment_aa.cigar is not None
            assert self.j_gene_alignment_aa.e_value is not None
            assert self.j_gene_alignment_aa.seq_identity is not None
            assert self.j_gene_alignment_aa.alignment_score is not None
            assert self.v_gene_alignment_aa.q_end is not None
            assert self.v_gene_alignment_aa.q_start is not None
            assert self.v_gene_alignment_aa.t_start is not None
            assert self.v_gene_alignment_aa.t_end is not None
            assert self.v_gene_alignment_aa.target_id is not None
            assert self.v_gene_alignment_aa.cigar is not None
            assert self.v_gene_alignment_aa.e_value is not None
            assert self.v_gene_alignment_aa.seq_identity is not None
            assert self.v_gene_alignment_aa.alignment_score is not None
            # offset c_aln by the end of j_gene_alignment_aa

            c_gene_alignment_aa = offset_alignments(self.j_gene_alignment_aa.q_end, c_aln)
            self.c_gene_alignment_aa = c_gene_alignment_aa
            self.rearrangement.c_call = c_gene_alignment_aa.target_id

            # alignment positions are 0-indexed, add 1 to start index to convert to 1-based
            self.rearrangement.c_sequence_start_aa = c_gene_alignment_aa.q_start + 1
            self.rearrangement.c_sequence_end_aa = c_gene_alignment_aa.q_end

            self.rearrangement.c_germline_start_aa = c_gene_alignment_aa.t_start + 1
            self.rearrangement.c_germline_end_aa = c_gene_alignment_aa.t_end

            c_sequence_segment = self.sequence_aa[c_gene_alignment_aa.q_start : c_gene_alignment_aa.q_end]
            c_germline_segment = c_aln.t_seq[c_aln.t_start : c_aln.t_end]

            c_sequence_alignment, c_germline_alignment = align_sequences(
                c_sequence_segment, c_germline_segment, unfold_cigar(c_aln.cigar)
            )

            self.rearrangement.c_sequence_alignment_aa = c_sequence_alignment
            self.rearrangement.c_germline_alignment_aa = c_germline_alignment

            self.rearrangement.c_score_aa = c_gene_alignment_aa.alignment_score
            self.rearrangement.c_cigar_aa = c_gene_alignment_aa.cigar
            self.rearrangement.c_support_aa = c_gene_alignment_aa.e_value
            self.rearrangement.c_identity_aa = c_gene_alignment_aa.seq_identity

        return self

    def _extract_region_aa(self, start: int, end: int):
        if start is None or end is None or start == -1 or end == -1:
            return None

        if self.rearrangement.sequence_aa:
            return self.rearrangement.sequence_aa[start - 1 : end]
        return None

    def with_aa_region_offsets(self, offsets: RegionOffsetsAA):
        self.aa_offsets = offsets

        self.rearrangement.fwr1_start_aa = offsets.fwr1_start_aa
        self.rearrangement.fwr1_end_aa = offsets.fwr1_end_aa
        self.rearrangement.cdr1_start_aa = offsets.cdr1_start_aa
        self.rearrangement.cdr1_end_aa = offsets.cdr1_end_aa
        self.rearrangement.fwr2_start_aa = offsets.fwr2_start_aa
        self.rearrangement.fwr2_end_aa = offsets.fwr2_end_aa
        self.rearrangement.cdr2_start_aa = offsets.cdr2_start_aa
        self.rearrangement.cdr2_end_aa = offsets.cdr2_end_aa
        self.rearrangement.fwr3_start_aa = offsets.fwr3_start_aa
        self.rearrangement.fwr3_end_aa = offsets.fwr3_end_aa
        self.rearrangement.cdr3_start_aa = offsets.cdr3_start_aa
        self.rearrangement.cdr3_end_aa = offsets.cdr3_end_aa
        self.rearrangement.fwr4_start_aa = offsets.fwr4_start_aa
        self.rearrangement.fwr4_end_aa = offsets.fwr4_end_aa

        self.rearrangement.fwr1_aa = self._extract_region_aa(offsets.fwr1_start_aa, offsets.fwr1_end_aa)
        self.rearrangement.cdr1_aa = self._extract_region_aa(offsets.cdr1_start_aa, offsets.cdr1_end_aa)
        self.rearrangement.fwr2_aa = self._extract_region_aa(offsets.fwr2_start_aa, offsets.fwr2_end_aa)
        self.rearrangement.cdr2_aa = self._extract_region_aa(offsets.cdr2_start_aa, offsets.cdr2_end_aa)
        self.rearrangement.fwr3_aa = self._extract_region_aa(offsets.fwr3_start_aa, offsets.fwr3_end_aa)
        self.rearrangement.cdr3_aa = self._extract_region_aa(offsets.cdr3_start_aa, offsets.cdr3_end_aa)
        self.rearrangement.fwr4_aa = self._extract_region_aa(offsets.fwr4_start_aa, offsets.fwr4_end_aa)

        return self

    def with_aa_scheme_alignment(self, alignment: SchemeAlignment):
        self.rearrangement.sequence_aa_scheme_cigar = fold_cigar(
            AlignmentString(alignment.alignment_str.replace("S", "").replace("N", ""))
        )
        return self

    def with_scheme_residue_mapping(self, mapping: dict[str, str]):
        self.rearrangement.scheme_residue_mapping = mapping
        return self

    def with_positional_scheme_mapping(self, mapping: dict[int, str]):
        self.rearrangement.positional_scheme_mapping = mapping
        return self

    def with_scheme_alignment_exc(self, exc: str):
        self.rearrangement.exc = exc
        return self

    def build(self):  # pylint: disable=too-many-statements
        if self.v_gene_alignment_aa is None:  # or self.rearrangement.sequence_alignment_aa is None:
            return self.rearrangement

        alignment_end = self.j_gene_alignment_aa.q_end if self.j_gene_alignment_aa else self.v_gene_alignment_aa.q_end
        sequence_alignment_aa_segment = self.sequence_aa[self.v_gene_alignment_aa.q_start : alignment_end]

        v_germline_segment = self.v_gene_alignment_aa.t_seq[
            self.v_gene_alignment_aa.t_start : self.v_gene_alignment_aa.t_end
        ]
        j_germline_segment = (
            self.j_gene_alignment_aa.t_seq[self.j_gene_alignment_aa.t_start : self.j_gene_alignment_aa.t_end]
            if self.j_gene_alignment_aa
            else ""
        )

        germline_segment = v_germline_segment + j_germline_segment

        v_aln_str = unfold_cigar(self.v_gene_alignment_aa.cigar)
        j_aln_str = unfold_cigar(self.j_gene_alignment_aa.cigar) if self.j_gene_alignment_aa else ""

        vj_junction_len = (
            self.j_gene_alignment_aa.q_start - self.v_gene_alignment_aa.q_end if self.j_gene_alignment_aa else 0
        )
        vj_junction_alignment_str = vj_junction_len * "I"

        full_alignment_str = v_aln_str + vj_junction_alignment_str + j_aln_str

        sequence_alignment_aa, germline_alignment_aa = align_sequences(
            sequence_alignment_aa_segment, germline_segment, full_alignment_str
        )
        self.rearrangement.sequence_alignment_aa = sequence_alignment_aa
        self.rearrangement.germline_alignment_aa = germline_alignment_aa

        if self.rearrangement.cdr3_aa is not None:
            jun_start_aa = self.aa_offsets.cdr3_start_aa - 1
            jun_end_aa = self.aa_offsets.cdr3_end_aa + 1

            self.rearrangement.junction_aa = self._extract_region_aa(jun_start_aa, jun_end_aa)
            self.rearrangement.junction_aa_length = len(self.rearrangement.junction_aa or "")

        # annotate
        self.rearrangement.stop_codon = "*" in self.rearrangement.sequence_alignment_aa
        if self.aa_offsets:
            self.rearrangement.complete_vdj = (
                self.aa_offsets.fwr1_start_aa
                and self.aa_offsets.fwr1_start_aa != -1
                and self.aa_offsets.fwr4_end_aa
                and self.aa_offsets.fwr4_end_aa != -1
            )

        self.rearrangement.productive = bool(
            not self.rearrangement.stop_codon and self.rearrangement.v_call and self.rearrangement.j_call
        )

        self.rearrangement.additional_validation_flags = validate_airr_entry_aa(self.rearrangement, scheme=self.scheme)

        return self.rearrangement


class SegmentedAirrBuilderAA(AirrBuilderAA):
    def __init__(self, sequence_header: str, sequence: str, scheme: Scheme, query_sequence: str):
        super().__init__(sequence_header, sequence, scheme)

        self.rearrangement: SegmentedAirrRearrangementEntryAA = SegmentedAirrRearrangementEntryAA(
            sequence_header=sequence_header,
            sequence_aa=sequence,
            numbering_scheme=scheme.value,
            query_sequence=query_sequence,
        )

    def with_segment_start_end(self, segment_start: int, segment_end: int):
        self.rearrangement.segment_start = segment_start + 1
        self.rearrangement.segment_end = segment_end
        return self
