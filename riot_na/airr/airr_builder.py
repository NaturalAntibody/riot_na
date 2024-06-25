from math import ceil
from typing import Optional

from riot_na.airr.airr_validator import validate_airr_entry
from riot_na.alignment.alignment_utils import fold_cigar, has_frameshift
from riot_na.data.model import (
    AirrRearrangementEntryNT,
    AlignmentEntryNT,
    AlignmentString,
    RegionOffsetsAA,
    RegionOffsetsNT,
    Scheme,
    SchemeAlignment,
)


class AirrBuilder:  # pylint: disable=too-many-instance-attributes
    def __init__(self, sequence_header: str, sequence: str, scheme: Scheme):
        self.sequence_header = sequence_header
        self.sequence = sequence
        self.scheme = scheme

        self.rearrangement = AirrRearrangementEntryNT(
            sequence_header=sequence_header, sequence=sequence, numbering_scheme=scheme.value
        )

        self.v_gene_sequence: Optional[str] = None
        self.v_gene_alignment: Optional[AlignmentEntryNT] = None

        self.j_gene_alignment: Optional[AlignmentEntryNT] = None
        self.j_gene_sequence: Optional[str] = None

        self.c_gene_alignment: Optional[AlignmentEntryNT] = None
        self.c_gene_sequence: Optional[str] = None
        self.c_gene_name: Optional[str] = None

        self.d_gene_alignment: Optional[AlignmentEntryNT] = None
        self.d_gene_sequence: Optional[str] = None
        self.d_gene_name: Optional[str] = None

        self.scheme_alingment: Optional[SchemeAlignment] = None
        self.nt_offsets: Optional[RegionOffsetsNT] = None
        self.aa_offsets: Optional[RegionOffsetsAA] = None

    def with_v_gene_alignment(
        self,
        v_aln: AlignmentEntryNT,
    ):
        self.v_gene_sequence = v_aln.t_seq
        self.rearrangement.rev_comp = v_aln.rev_comp
        self.v_gene_alignment = v_aln

        if v_aln is None:
            return self

        self.rearrangement.v_call = v_aln.target_id
        self.rearrangement.locus = v_aln.locus.value

        # alignment positions are 0-indexed, add 1 to start index to convert to 1-based
        self.rearrangement.v_sequence_start = v_aln.q_start + 1
        self.rearrangement.v_sequence_end = v_aln.q_end

        self.rearrangement.v_germline_start = v_aln.t_start + 1
        self.rearrangement.v_germline_end = v_aln.t_end

        self.rearrangement.v_sequence_alignment = self.sequence[v_aln.q_start : v_aln.q_end]

        # 1 based, hence need to add 1
        self.rearrangement.v_alignment_start = 1
        self.rearrangement.v_alignment_end = v_aln.q_end - v_aln.q_start

        self.rearrangement.v_germline_alignment = v_aln.t_seq[v_aln.t_start : v_aln.t_end]

        self.rearrangement.v_score = v_aln.alignment_score
        self.rearrangement.v_cigar = v_aln.cigar
        self.rearrangement.v_support = v_aln.e_value
        self.rearrangement.v_identity = v_aln.seq_identity

        return self

    def with_j_gene_alignment(
        self,
        j_aln: AlignmentEntryNT,
    ):
        self.j_gene_sequence = j_aln.t_seq

        if j_aln is not None:
            assert self.v_gene_alignment is not None

            self.j_gene_alignment = j_aln

            self.rearrangement.j_call = j_aln.target_id
            # alignment positions are 0-indexed, add 1 to start index to convert to 1-based
            self.rearrangement.j_sequence_start = j_aln.q_start + 1
            self.rearrangement.j_sequence_end = j_aln.q_end

            self.rearrangement.j_germline_start = j_aln.t_start + 1
            self.rearrangement.j_germline_end = j_aln.t_end

            self.rearrangement.j_sequence_alignment = self.sequence[j_aln.q_start : j_aln.q_end]

            # 1 based, hence need to add 1
            assert self.rearrangement.v_sequence_start is not None
            self.rearrangement.j_alignment_start = j_aln.q_start + 1 - self.v_gene_alignment.q_start
            self.rearrangement.j_alignment_end = j_aln.q_end - self.v_gene_alignment.q_start

            self.rearrangement.j_germline_alignment = j_aln.t_seq[j_aln.t_start : j_aln.t_end]

            self.rearrangement.j_score = j_aln.alignment_score
            self.rearrangement.j_cigar = j_aln.cigar
            self.rearrangement.j_support = j_aln.e_value
            self.rearrangement.j_identity = j_aln.seq_identity

            assert j_aln.reading_frame is not None
            self.rearrangement.j_frame = (3 - ((j_aln.t_start - j_aln.reading_frame) % 3)) % 3

        return self

    def with_c_gene_alignment(
        self,
        c_aln: AlignmentEntryNT,
    ):
        self.c_gene_sequence = c_aln.t_seq
        self.c_gene_name = c_aln.target_id

        if self.j_gene_alignment is None:
            return self

        self.c_gene_alignment = c_aln

        self.rearrangement.c_call = c_aln.target_id
        # alignment positions are 0-indexed
        self.rearrangement.c_sequence_start = c_aln.q_start + 1
        self.rearrangement.c_sequence_end = c_aln.q_end

        self.rearrangement.c_germline_start = c_aln.t_start + 1
        self.rearrangement.c_germline_end = c_aln.t_end

        self.rearrangement.c_sequence_alignment = self.sequence[c_aln.q_start : c_aln.q_end]

        # 1 based, hence need to add 1
        assert self.v_gene_alignment is not None
        self.rearrangement.c_alignment_start = c_aln.q_start + 1 - self.v_gene_alignment.q_start
        self.rearrangement.c_alignment_end = c_aln.q_end - self.v_gene_alignment.q_start

        self.rearrangement.c_germline_alignment = c_aln.t_seq[c_aln.t_start : c_aln.t_end]

        self.rearrangement.c_score = c_aln.alignment_score
        self.rearrangement.c_cigar = c_aln.cigar
        self.rearrangement.c_support = c_aln.e_value
        self.rearrangement.c_identity = c_aln.seq_identity

        return self

    def with_d_gene_alignment(
        self,
        d_aln: AlignmentEntryNT,
    ):
        self.d_gene_sequence = d_aln.t_seq
        self.d_gene_name = d_aln.target_id

        if self.j_gene_alignment is None or self.v_gene_alignment is None:
            return self

        self.d_gene_alignment = d_aln

        self.rearrangement.d_call = d_aln.target_id
        # alignment positions are 0-indexed
        self.rearrangement.d_sequence_start = d_aln.q_start + 1
        self.rearrangement.d_sequence_end = d_aln.q_end

        self.rearrangement.d_germline_start = d_aln.t_start + 1
        self.rearrangement.d_germline_end = d_aln.t_end

        self.rearrangement.d_sequence_alignment = self.sequence[d_aln.q_start : d_aln.q_end]

        assert self.v_gene_alignment is not None
        self.rearrangement.d_alignment_start = d_aln.q_start + 1 - self.v_gene_alignment.q_start
        self.rearrangement.d_alignment_end = d_aln.q_end - self.v_gene_alignment.q_start

        self.rearrangement.d_germline_alignment = d_aln.t_seq[d_aln.t_start : d_aln.t_end]

        self.rearrangement.d_score = d_aln.alignment_score
        self.rearrangement.d_cigar = d_aln.cigar
        self.rearrangement.d_support = d_aln.e_value
        self.rearrangement.d_identity = d_aln.seq_identity

        return self

    def with_sequence_alignment_aa(self, aligned_sequence_aa: Optional[str]):
        self.rearrangement.sequence_alignment_aa = aligned_sequence_aa
        return self

    def with_sequence_alignment(self, aligned_sequence: str):
        self.rearrangement.sequence_alignment = aligned_sequence
        return self

    def with_reading_frame(self, reading_frame: int):
        self.rearrangement.v_frame = reading_frame
        return self

    def with_nt_region_offsets(self, offsets: RegionOffsetsNT):
        self.nt_offsets = offsets

        self.rearrangement.fwr1_start = offsets.fwr1_start
        self.rearrangement.fwr1_end = offsets.fwr1_end
        self.rearrangement.cdr1_start = offsets.cdr1_start
        self.rearrangement.cdr1_end = offsets.cdr1_end
        self.rearrangement.fwr2_start = offsets.fwr2_start
        self.rearrangement.fwr2_end = offsets.fwr2_end
        self.rearrangement.cdr2_start = offsets.cdr2_start
        self.rearrangement.cdr2_end = offsets.cdr2_end
        self.rearrangement.fwr3_start = offsets.fwr3_start
        self.rearrangement.fwr3_end = offsets.fwr3_end
        self.rearrangement.cdr3_start = offsets.cdr3_start
        self.rearrangement.cdr3_end = offsets.cdr3_end
        self.rearrangement.fwr4_start = offsets.fwr4_start
        self.rearrangement.fwr4_end = offsets.fwr4_end

        self.rearrangement.fwr1 = self._extract_region(offsets.fwr1_start, offsets.fwr1_end)
        self.rearrangement.cdr1 = self._extract_region(offsets.cdr1_start, offsets.cdr1_end)
        self.rearrangement.fwr2 = self._extract_region(offsets.fwr2_start, offsets.fwr2_end)
        self.rearrangement.cdr2 = self._extract_region(offsets.cdr2_start, offsets.cdr2_end)
        self.rearrangement.fwr3 = self._extract_region(offsets.fwr3_start, offsets.fwr3_end)
        self.rearrangement.cdr3 = self._extract_region(offsets.cdr3_start, offsets.cdr3_end)
        self.rearrangement.fwr4 = self._extract_region(offsets.fwr4_start, offsets.fwr4_end)

        return self

    def with_aa_region_offsets(self, offsets: RegionOffsetsAA):
        self.aa_offsets = offsets

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

    def _extract_region(self, start: int, end: int):
        if start is None or end is None or start == -1 or end == -1:
            return None

        assert self.v_gene_alignment is not None
        start = start - self.v_gene_alignment.q_start
        end = end - self.v_gene_alignment.q_start

        if self.rearrangement.sequence_alignment:
            # subtract 1 because of 1 based indexing
            return self.rearrangement.sequence_alignment[start - 1 : end]
        return None

    def _extract_region_aa(self, start: int, end: int):
        if start is None or end is None or start == -1 or end == -1:
            return None

        if self.rearrangement.sequence_alignment_aa:
            return self.rearrangement.sequence_alignment_aa[start - 1 : end]
        return None

    def build(self):  # pylint: disable=too-many-statements
        if self.v_gene_alignment is None:
            return self.rearrangement

        alignment_end = self.j_gene_alignment.q_end if self.j_gene_alignment else self.v_gene_alignment.q_end
        self.rearrangement.sequence_alignment = self.sequence[self.v_gene_alignment.q_start : alignment_end]

        aln_len = len(self.rearrangement.sequence_alignment)
        v_germline_alignment = self.rearrangement.v_germline_alignment
        j_germline_alignment = self.rearrangement.j_germline_alignment or ""

        self.rearrangement.germline_alignment = (
            v_germline_alignment
            + ("N" * (aln_len - len(v_germline_alignment) - len(j_germline_alignment)))
            + j_germline_alignment
        )

        if self.rearrangement.sequence_alignment_aa:
            alignment_offset = self.v_gene_alignment.q_start + self.rearrangement.v_frame
            v_start_aa = 0
            v_end_aa = int((self.v_gene_alignment.q_end - alignment_offset) / 3)

            self.rearrangement.v_sequence_alignment_aa = self.rearrangement.sequence_alignment_aa[v_start_aa:v_end_aa]

            if self.j_gene_alignment is not None:
                j_start_aa = ceil((self.j_gene_alignment.q_start - alignment_offset) / 3)
                j_end_aa = int(self.j_gene_alignment.q_end - alignment_offset / 3)

                self.rearrangement.j_sequence_alignment_aa = self.rearrangement.sequence_alignment_aa[
                    j_start_aa:j_end_aa
                ]

            if self.rearrangement.cdr3 is not None:
                self.rearrangement.junction = self._extract_region(
                    self.rearrangement.cdr3_start - 3, self.rearrangement.cdr3_end + 3
                )
                self.rearrangement.junction_length = len(self.rearrangement.junction or "")

            if self.rearrangement.cdr3_aa is not None:
                jun_start_aa = self.aa_offsets.cdr3_start_aa - 1
                jun_end_aa = self.aa_offsets.cdr3_end_aa + 1

                self.rearrangement.junction_aa = self._extract_region_aa(jun_start_aa, jun_end_aa)
                self.rearrangement.junction_aa_length = len(self.rearrangement.junction_aa or "")

            # annotate
            self.rearrangement.stop_codon = "*" in self.rearrangement.sequence_alignment_aa
        self.rearrangement.complete_vdj = (
            self.rearrangement.fwr1_start
            and self.rearrangement.fwr1_start != -1
            and self.rearrangement.fwr4_end
            and self.rearrangement.fwr4_end != -1
        )

        # airr defines frameshift as sum of all insertions and deletions on gene alignment
        self.rearrangement.v_frameshift = self.rearrangement.v_cigar and has_frameshift(self.rearrangement.v_cigar)
        self.rearrangement.j_frameshift = self.rearrangement.j_cigar and has_frameshift(self.rearrangement.j_cigar)

        if self.rearrangement.j_frame is not None:
            j_frame_start = self.rearrangement.j_alignment_start - 1 + self.rearrangement.j_frame
            self.rearrangement.vj_in_frame = (j_frame_start % 3) == self.rearrangement.v_frame

            productive = not self.rearrangement.stop_codon and self.rearrangement.vj_in_frame
            self.rearrangement.productive = productive
        else:
            self.rearrangement.vj_in_frame = None
            self.rearrangement.productive = False

        self.rearrangement.additional_validation_flags = validate_airr_entry(self.rearrangement, scheme=self.scheme)

        return self.rearrangement
