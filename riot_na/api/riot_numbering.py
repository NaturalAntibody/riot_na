# pylint: disable=too-many-branches, too-many-statements

from pathlib import Path
from typing import Optional

from riot_na.airr.airr_builder import AirrBuilder
from riot_na.airr.airr_builder_aa import AirrBuilderAA
from riot_na.alignment.aa_gene_alignments import (
    VJAlignerAA,
    VJAlignmentTranslatorAA,
    create_vj_aligner_aa,
    create_vj_alignment_translator_aa,
)
from riot_na.alignment.nt_gene_alignments import VDJCAlignerNT, create_vdjc_aligner_nt
from riot_na.config import GENE_DB_DIR
from riot_na.data.model import (
    AirrRearrangementEntryAA,
    AirrRearrangementEntryNT,
    Organism,
    Scheme,
)
from riot_na.schemes.region_offsets import (
    infer_aa_region_offsets,
    infer_nt_region_offsets,
    infer_region_offsets,
)
from riot_na.schemes.scheme_alignment import (
    SchemeAligner,
    fix_imgt_cdrs_numbering,
    get_positional_scheme_mapping,
    get_scheme_residue_mapping,
    produce_numbering,
    sort_imgt_numbering,
)


class RiotNumberingNT:
    def __init__(
        self,
        vdjc_aligner_nt: VDJCAlignerNT,
        vj_alignment_translator_aa: VJAlignmentTranslatorAA,
        scheme_aligner: SchemeAligner,
    ):
        self.vdjc_aligner = vdjc_aligner_nt
        self.vj_alignment_translator_aa = vj_alignment_translator_aa
        self.scheme_aligner = scheme_aligner

    def run_on_sequence(
        self,
        header: str,
        query_sequence: str,
        scheme: Scheme = Scheme.IMGT,
    ) -> AirrRearrangementEntryNT:  # pylint: disable=too-many-statements
        """
        Align and number nucleotide sequence.

        Parameters:
        -----------
        - header : Fasta sequence header, usually this is simply sequence ID.
        - query_sequence : Nucleotide sequence to number.
        - scheme : Which numbering scheme to use (default IMGT).

        Returns:
        -----------
        AirrRearrangementEntryNT object with alignment results (Extended AIRR format)
        """
        query_sequence = query_sequence.upper()
        translated_alignments = None
        nt_offsets = None
        aa_offsets = None
        sch_alignment = None
        positional_scheme_mapping = None
        numbering = None
        scheme_alignment_exc = None

        # produce alignments for v,d,j,c
        alignments = self.vdjc_aligner.produce_nt_alignments(query_sequence)

        # Query might have been reverse complemented during V alignment.
        if alignments.v:
            query_sequence = alignments.v.q_seq
        try:
            translated_alignments = self.vj_alignment_translator_aa.translate_nt_alignments(alignments)

            if translated_alignments:
                sch_alignment = self.scheme_aligner.align_to_scheme(translated_alignments.aa_alignments, scheme)

                # calculate aa offsets based on the numbering (iterate over cigar to produce positional-> imgt mapping)
                assert sch_alignment
                assert alignments.v
                locus = alignments.v.locus

                region_offsets = infer_region_offsets(sch_alignment, scheme, locus)
                aa_offsets = infer_aa_region_offsets(region_offsets)

                alignment_start = alignments.v.q_start

                nt_offsets = infer_nt_region_offsets(
                    region_offsets, alignment_start, translated_alignments.reading_frame
                )
                numbering = produce_numbering(
                    translated_alignments.translated_query, sch_alignment.alignment_str, sch_alignment.q_start
                )

                if scheme == Scheme.IMGT:
                    numbering = fix_imgt_cdrs_numbering(numbering)
                    numbering = sort_imgt_numbering(numbering)

                scheme_residue_mapping = get_scheme_residue_mapping(numbering)
                positional_scheme_mapping = get_positional_scheme_mapping(scheme_residue_mapping, sch_alignment.q_start)
        except AssertionError as err:
            scheme_alignment_exc = str(err)

        # produce airr result
        # this is separated on purpose to enable other output formats
        airr_builder = AirrBuilder(header, query_sequence, scheme)

        if alignments.v:
            airr_builder.with_v_gene_alignment(alignments.v)

        if translated_alignments:
            airr_builder.with_aa_alignments(translated_alignments)

        if alignments.j:
            airr_builder.with_j_gene_alignment(alignments.j)

        if alignments.d:
            airr_builder.with_d_gene_alignment(alignments.d)

        if alignments.c:
            airr_builder.with_c_gene_alignment(alignments.c)

        if sch_alignment:
            assert nt_offsets
            assert aa_offsets
            airr_builder.with_nt_region_offsets(nt_offsets)
            airr_builder.with_aa_region_offsets(aa_offsets)
            airr_builder.with_aa_scheme_alignment(sch_alignment)

        if numbering:
            (
                airr_builder.with_scheme_residue_mapping(scheme_residue_mapping).with_positional_scheme_mapping(
                    positional_scheme_mapping
                )
            )

        if scheme_alignment_exc:
            airr_builder.with_scheme_alignment_exc(scheme_alignment_exc)

        result = airr_builder.build()
        return result


class RiotNumberingAA:
    def __init__(self, vj_aligner_aa: VJAlignerAA, scheme_aligner: SchemeAligner):
        self.vj_aligner_aa = vj_aligner_aa
        self.scheme_aligner = scheme_aligner

    def run_on_sequence(
        self,
        header: str,
        query_sequence: str,
        scheme: Scheme = Scheme.IMGT,
        extend_alignment: bool = True,
    ) -> AirrRearrangementEntryAA:  # pylint: disable=too-many-statements
        """
        Align and number amino acid sequence.

        Parameters:
        -----------
        - header : Fasta sequence header, usually this is simply sequence ID.
        - query_sequence : Nucleotide sequence to number.
        - scheme : Which numbering scheme to use (default IMGT).
        - extend_alignment: RIOT uses Striped Smith-Waterman local alignemtn algorithm.
        Due to that, if query sequence has too many differences (due to mutations or other reasons)
        to the germline sequence at 5' or / and 3' ends the result may not show that part of the gene.
        This option appends missing amino acids, so the result covers full gene sequence. Beware that
        those query fragments are not truly a part of alignemnt and are shown only for convenience!

        Returns:
        -----------
        AirrRearrangementEntryAA object with alignment results (Extended AIRR format)
        """
        query_sequence = query_sequence.upper()
        sch_alignment = None
        aa_offsets = None
        aligned_sequence = None
        positional_scheme_mapping = None
        numbering = None
        scheme_alignment_exc = None

        # produce alignments for v,d,j,c
        alignments = self.vj_aligner_aa.produce_aa_alignments(query_sequence)

        if alignments.v:
            try:
                sch_alignment = self.scheme_aligner.align_to_scheme(
                    alignments, scheme, extend_alignment=extend_alignment
                )

                # calculate aa offsets based on the numbering (iterate over cigar to produce positional-> imgt mapping)
                assert sch_alignment
                locus = alignments.v.locus

                region_offsets = infer_region_offsets(sch_alignment, scheme, locus)
                aa_offsets = infer_aa_region_offsets(region_offsets)

                alignment_start = sch_alignment.q_start
                alignment_end = sch_alignment.q_end

                aligned_sequence = query_sequence[alignment_start:alignment_end]

                numbering = produce_numbering(aligned_sequence, sch_alignment.alignment_str)

                if scheme == Scheme.IMGT:
                    numbering = fix_imgt_cdrs_numbering(numbering)
                    numbering = sort_imgt_numbering(numbering)

                scheme_residue_mapping = get_scheme_residue_mapping(numbering)
                positional_scheme_mapping = get_positional_scheme_mapping(scheme_residue_mapping, sch_alignment.q_start)
            except AssertionError as err:
                scheme_alignment_exc = str(err)

        # produce airr result
        # this is separated on purpose to enable other output formats
        airr_builder = AirrBuilderAA(header, query_sequence, scheme)

        if alignments.v:
            airr_builder.with_v_gene_alignment_aa(alignments.v)

        if alignments.j:
            airr_builder.with_j_gene_alignment_aa(alignments.j)

        if sch_alignment:
            assert aa_offsets
            airr_builder.with_aa_region_offsets(aa_offsets)
            airr_builder.with_aa_scheme_alignment(sch_alignment)

        if numbering:
            (
                airr_builder.with_scheme_residue_mapping(scheme_residue_mapping).with_positional_scheme_mapping(
                    positional_scheme_mapping
                )
            )

        if scheme_alignment_exc:
            airr_builder.with_scheme_alignment_exc(scheme_alignment_exc)

        result = airr_builder.build()
        return result


def create_riot_nt(
    allowed_species: Optional[list[Organism]] = None,
    db_dir: Path = GENE_DB_DIR,
) -> RiotNumberingNT:
    """
    Factory function for creating RiotNumberingNT object.

    Parameters:
    -----------
    - allowed_species : Limit gene database to specified Organism list (default use all species - mouse and human).
    - allowed_schemes : Limit numbering capabilities to the list of specified schemes (default use all schemes).
    - db_dir : Path to gene and scheme mappings database directory (default use embeded database).

    Returns:
    -----------
    RiotNumberingNT object used for performing alignment and numbering.
    """
    vdjc_aligner_nt = create_vdjc_aligner_nt(allowed_species=allowed_species, db_dir=db_dir)
    vj_alignment_translator_aa = create_vj_alignment_translator_aa(allowed_species=allowed_species, db_dir=db_dir)
    scheme_aligner = SchemeAligner(allowed_species=allowed_species, db_dir=db_dir)
    return RiotNumberingNT(vdjc_aligner_nt, vj_alignment_translator_aa, scheme_aligner)


def create_riot_aa(
    allowed_species: Optional[list[Organism]] = None,
    db_dir: Path = GENE_DB_DIR,
) -> RiotNumberingAA:
    """
    Factory function for creating RiotNumberingNT object.

    Parameters:
    -----------
    - allowed_species : Limit gene database to specified Organism list (default use all species - mouse and human).
    - allowed_schemes : Limit numbering capabilities to the list of specified schemes (default use all schemes).
    - db_dir : Path to gene and scheme mappings database directory (default use embeded database).

    Returns:
    -----------
    RiotNumberingAA object used for performing alignment and numbering.
    """
    vj_aligner_aa = create_vj_aligner_aa(
        allowed_species=allowed_species,
        db_dir=db_dir,
    )
    scheme_aligner = SchemeAligner(allowed_species=allowed_species, db_dir=db_dir)
    return RiotNumberingAA(vj_aligner_aa, scheme_aligner)


if __name__ == "__main__":
    import json
    from dataclasses import asdict

    vdj_alnr = create_vdjc_aligner_nt()
    vdj_aa_translator = create_vj_alignment_translator_aa()
    scheme_alnr = SchemeAligner()

    # riot_numbering = RiotNumberingNT(vdj_alnr, vdj_aa_translator, scheme_alnr)
    # # given query
    # QUERY = "CAGGTGCAGCTACAGCAGTGGGGCGCAGGACTGTTGAAGCCTTCGGAGACCCTGTCCCTCACCTGCGCTGTCTATGGTGGGTCCTTCAGTGGTTACTACTGGAGCTGGATCCGCCAGCCCCCAGGGAAGGGGCTGGAGTGGATTGGGGAAATCAATCATAGTGGAAGCACCAACTACAACCCGTCCCTCAAGAGTCGAGTCACCATATCAGTAGACACGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCGGACACGGCTGTGTATTACTGTGCGAGAGAAAGGCGAGGTTGTAGTAGTACCAGCTGCTTTTGGATTACTCAACGTGGATACAGCTATGGAGGCTCCTACTACTACTACTACATGGACGTCTGGGGCAAAGGGACCACGGTCACCGTCTCCTCA"
    # sample_result = riot_numbering.run_on_sequence("header", QUERY, Scheme.IMGT)
    # print(json.dumps(asdict(sample_result), indent=4))

    vdj_aa_alnr = create_vj_aligner_aa()
    aa_numbering = RiotNumberingAA(vdj_aa_alnr, scheme_alnr)

    AA_QUERY = "QLQLQESGPGLVKPSETLSLTCTVSGDSISSGDYYWGWIRQPPGKGLEWIGHIYYSGATYYNPSLENRVTISVDTSKNQFSLKLSSVTAADTAVYYCTRDDSSNWRSRGQGTLVTVSS"

    aa_sample_result = aa_numbering.run_on_sequence("header", AA_QUERY, Scheme.IMGT)
    print(json.dumps(asdict(aa_sample_result), indent=4))
