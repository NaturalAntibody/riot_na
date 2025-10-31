# pylint: disable=too-many-branches, too-many-statements

from pathlib import Path
from typing import Optional

from cachetools import cached

from riot_na.airr.airr_builder import AirrBuilder, SegmentedAirrBuilder
from riot_na.airr.airr_builder_aa import AirrBuilderAA, SegmentedAirrBuilderAA
from riot_na.alignment.aa_gene_alignments import (
    VJCAlignerAA,
    VJCAlignmentTranslatorAA,
    create_vjc_aligner_aa,
    create_vjc_alignment_translator_aa,
)
from riot_na.alignment.nt_gene_alignments import VDJCAlignerNT, create_vdjc_aligner_nt
from riot_na.config import GENE_DB_DIR
from riot_na.data.model import (
    AirrRearrangementEntryAA,
    AirrRearrangementEntryNT,
    ChainType,
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
        vjc_alignment_translator_aa: VJCAlignmentTranslatorAA,
        scheme_aligner: SchemeAligner,
        return_all_domains: bool = False,
    ):
        self.vdjc_aligner = vdjc_aligner_nt
        self.vjc_alignment_translator_aa = vjc_alignment_translator_aa
        self.scheme_aligner = scheme_aligner
        self.return_all_domains = return_all_domains

    def run_on_sequence(
        self,
        header: str,
        query_sequence: str,
        scheme: Scheme = Scheme.IMGT,
        extend_alignment: bool = True,
    ) -> AirrRearrangementEntryNT | list[AirrRearrangementEntryNT]:  # pylint: disable=too-many-statements
        """
        Align and number nucleotide sequence.

        Parameters:
        -----------
        - header : Fasta sequence header, usually this is simply sequence ID.
        - query_sequence : Nucleotide sequence to number.
        - scheme : Which numbering scheme to use (default IMGT).
        - return_all_domains: If True, return all domains of multiple domain proteins.
        - extend_alignment: Extend alignments to cover full sequence.

        Returns:
        -----------
        AirrRearrangementEntryNT object with alignment results (Extended AIRR format)
        """
        query_sequence = query_sequence.upper().replace("U", "T")
        sequence = ""
        segment_start = -1
        segment_end = -1
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
        results = []

        for alignment in alignments:
            if alignment.v:
                sequence = alignment.v.q_seq
                segment_start = alignment.segment_start or 0
                segment_end = alignment.segment_end or len(sequence)
                try:
                    translated_alignments = self.vjc_alignment_translator_aa.translate_nt_alignments(
                        alignment, extend_alignment=extend_alignment
                    )

                    if translated_alignments:
                        sch_alignment = self.scheme_aligner.align_to_scheme(
                            translated_alignments.aa_alignments, scheme, extend_alignment=False
                        )

                        # calculate aa offsets based on the numbering (iterate over cigar to produce positional-> imgt mapping)
                        assert sch_alignment
                        assert alignment.v
                        chain_type = ChainType.from_locus(alignment.v.locus)

                        region_offsets = infer_region_offsets(sch_alignment, scheme, chain_type)
                        aa_offsets = infer_aa_region_offsets(region_offsets)

                        alignment_start = alignment.v.q_start

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
                        positional_scheme_mapping = get_positional_scheme_mapping(
                            scheme_residue_mapping, sch_alignment.q_start
                        )
                except AssertionError as err:
                    scheme_alignment_exc = str(err)

            # produce airr result
            # this is separated on purpose to enable other output formats
            airr_builder = AirrBuilder(header, sequence or query_sequence, scheme)
            if self.return_all_domains:
                airr_builder = SegmentedAirrBuilder(header, sequence, scheme, query_sequence)
                airr_builder.with_segment_start_end(segment_start, segment_end)

            if alignment.v:
                airr_builder.with_v_gene_alignment(alignment.v)

            if translated_alignments:
                airr_builder.with_aa_alignments(translated_alignments)

            if alignment.j:
                airr_builder.with_j_gene_alignment(alignment.j)

            if alignment.d:
                airr_builder.with_d_gene_alignment(alignment.d)

            if alignment.c:
                airr_builder.with_c_gene_alignment(alignment.c)

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
            results.append(result)

        if not results:
            results = [
                (
                    AirrBuilder(header, sequence, scheme).build()
                    if not self.return_all_domains
                    else SegmentedAirrBuilder(header, sequence, scheme, query_sequence).build()
                )
            ]

        if self.return_all_domains:
            return results

        return sorted(results)[-1]


class RiotNumberingAA:
    def __init__(self, vjc_aligner_aa: VJCAlignerAA, scheme_aligner: SchemeAligner, return_all_domains: bool = False):
        self.vjc_aligner_aa = vjc_aligner_aa
        self.scheme_aligner = scheme_aligner
        self.return_all_domains = return_all_domains

    def run_on_sequence(
        self,
        header: str,
        query_sequence: str,
        scheme: Scheme = Scheme.IMGT,
        extend_alignment: bool = True,
    ) -> AirrRearrangementEntryAA | list[AirrRearrangementEntryAA]:  # pylint: disable=too-many-statements
        """
        Align and number amino acid sequence.

        Parameters:
        -----------
        - header : Fasta sequence header, usually this is simply sequence ID.
        - query_sequence : Amino acid sequence to number.
        - scheme : Which numbering scheme to use (default IMGT).
        - extend_alignment: RIOT uses Striped Smith-Waterman local alignemtn algorithm.
        Due to that, if query sequence has too many differences (due to mutations or other reasons)
        to the germline sequence at 5' or / and 3' ends the result may not show that part of the gene.
        This option appends missing amino acids, so the result covers full gene sequence. Beware that
        those query fragments are not truly a part of alignemnt and are shown only for convenience!
        - return_all_domains: If True, return all domains of multiple domain proteins.

        Returns:
        -----------
        AirrRearrangementEntryAA object with alignment results (Extended AIRR format)
        """
        query_sequence = query_sequence.upper()
        sequence_aa = query_sequence
        sch_alignment = None
        aa_offsets = None
        positional_scheme_mapping = None
        numbering = None
        scheme_alignment_exc = None
        segment_start = -1
        segment_end = -1

        # produce alignments for v,j,c
        alignments = self.vjc_aligner_aa.produce_aa_alignments(query_sequence)

        results = []

        for alignment in alignments:
            if alignment.v:
                segment_start = alignment.segment_start or 0
                segment_end = alignment.segment_end or len(alignment.v.q_seq)

                sequence_aa = alignment.v.q_seq
                try:
                    sch_alignment = self.scheme_aligner.align_to_scheme(
                        alignment, scheme, extend_alignment=extend_alignment
                    )

                    # calculate aa offsets based on the numbering (iterate over cigar to produce positional-> imgt mapping)
                    assert sch_alignment
                    chain_type = ChainType.from_locus(alignment.v.locus)

                    region_offsets = infer_region_offsets(sch_alignment, scheme, chain_type)
                    aa_offsets = infer_aa_region_offsets(region_offsets)

                    numbering = produce_numbering(alignment.v.q_seq, sch_alignment.alignment_str, sch_alignment.q_start)

                    if scheme == Scheme.IMGT:
                        numbering = fix_imgt_cdrs_numbering(numbering)
                        numbering = sort_imgt_numbering(numbering)

                    scheme_residue_mapping = get_scheme_residue_mapping(numbering)
                    positional_scheme_mapping = get_positional_scheme_mapping(
                        scheme_residue_mapping, sch_alignment.q_start
                    )
                except AssertionError as err:
                    scheme_alignment_exc = str(err)

            # produce airr result
            # this is separated on purpose to enable other output formats
            airr_builder = AirrBuilderAA(header, query_sequence, scheme)
            if self.return_all_domains:
                airr_builder = SegmentedAirrBuilderAA(header, sequence_aa, scheme, query_sequence)
                airr_builder.with_segment_start_end(segment_start, segment_end)

            if alignment.v:
                airr_builder.with_v_gene_alignment_aa(alignment.v)

            if alignment.j:
                airr_builder.with_j_gene_alignment_aa(alignment.j)

            if alignment.c:
                airr_builder.with_c_gene_alignment_aa(alignment.c)

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
            results.append(result)

        if not results:
            results = [
                (
                    AirrBuilderAA(header, sequence_aa, scheme).build()
                    if not self.return_all_domains
                    else SegmentedAirrBuilderAA(header, sequence_aa, scheme, query_sequence).build()
                )
            ]

        if self.return_all_domains:
            return results

        return sorted(results)[-1]


def create_riot_nt(
    allowed_species: Optional[list[Organism]] = None,
    return_all_domains: bool = False,
    db_dir: Path = GENE_DB_DIR,
) -> RiotNumberingNT:
    """
    Factory function for creating RiotNumberingNT object.

    Parameters:
    -----------
    - allowed_species : Limit gene database to specified Organism list (default use all species - mouse and human).
    - return_all_domains: If True, return all domains of multiple domain proteins.
    - db_dir : Path to gene and scheme mappings database directory (default use embeded database).

    Returns:
    -----------
    RiotNumberingNT object used for performing alignment and numbering.
    """
    vdjc_aligner_nt = create_vdjc_aligner_nt(
        allowed_species=allowed_species, db_dir=db_dir, use_segment_aligner=return_all_domains
    )
    vjc_alignment_translator_aa = create_vjc_alignment_translator_aa(allowed_species=allowed_species, db_dir=db_dir)
    scheme_aligner = SchemeAligner(allowed_species=allowed_species, db_dir=db_dir)
    return RiotNumberingNT(vdjc_aligner_nt, vjc_alignment_translator_aa, scheme_aligner, return_all_domains)


def create_riot_aa(
    allowed_species: Optional[list[Organism]] = None,
    return_all_domains: bool = False,
    db_dir: Path = GENE_DB_DIR,
) -> RiotNumberingAA:
    """
    Factory function for creating RiotNumberingAA object.

    Parameters:
    -----------
    - allowed_species : Limit gene database to specified Organism list (default use all species - mouse and human).
    - return_all_domains: If True, return all domains of multiple domain proteins.
    - db_dir : Path to gene and scheme mappings database directory (default use embeded database).

    Returns:
    -----------
    RiotNumberingAA object used for performing alignment and numbering.
    """
    vjc_aligner_aa = create_vjc_aligner_aa(
        allowed_species=allowed_species,
        db_dir=db_dir,
        use_segment_aligner=return_all_domains,
    )
    scheme_aligner = SchemeAligner(allowed_species=allowed_species, db_dir=db_dir)
    return RiotNumberingAA(vjc_aligner_aa, scheme_aligner, return_all_domains)


@cached({})
def get_or_create_riot_nt(**riot_kwargs) -> RiotNumberingNT:
    return create_riot_nt(**riot_kwargs)


@cached({})
def get_or_create_riot_aa(**riot_kwargs) -> RiotNumberingAA:
    return create_riot_aa(**riot_kwargs)


if __name__ == "__main__":
    import json
    from dataclasses import asdict

    nt_numbering = create_riot_nt(return_all_domains=True)
    NT_QUERY = "caagugcaauuggucgagagcgggggcggggucgugcaaccgggcaggagccugcgcuugagcugugccgcgagccaguucacguucgggagcuauggcaugcacugggucaggcaaauacccgggaaagggcucgaguggguggccacgaucucguacgauggcacgaaaaaguaccaugccgauagcgugugggaucgcuuuaucauaagcagggacaauagcaagaacacgcucuuucuccagaugaauagcuugcgcccggaagauaccgcgcucuauuucugugucaaggaccagcgccaggacgagugcgaggaaugguggucggacuauuacgauuucgggcgcaggcucccgugcaggaaaucgcgcgggcuagcgggcauauucgaugucuggggccauggcacgauggugaccgugagcucg"
    # NT_QUERY = "CAGGCCCACCTGGAGCAAAGCGGCTCCGGGGTGAAGAAACCCGGAGCTTCTGTCAGAGTTAGCTGCTGGTCCTCTGAAGACATCTTCGAGCGGACCGAACTCATTCATTGGGTGCGCCAGGCCCCTGGCCAGGGGCTGGAGTGGATCGGATGGGTGAAGACAGTCACGGGCGCGGTGAACTTTGGGAGCCTTGATTTCAGACACCGGATTTCCCTGACCCGCGACAGAGATCTCTCTACAGCTTACATGGACATCCGGGGACTGACCCAGGATGACACAGCCACCTATTTTTGTGCCCGCCAAAAATTCGCTAGCAGATACTCCGGCGATCAGGGCAGCTATTTTGACCTTTGGGGGCGGGGAACACTGATTGTTGTGTCTTCC"
    # NT_QUERY = (
    #         # VH domain (nucleotide)
    #         "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTCCAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGCAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGCAAGGGGCTGGAGTGGGTGGCAGTTATATCATATGATGGAAGTAATAAATACTATGCAGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGAGATCGGGGGAGGAACGGTTATGATGCTTTTGATATCTGGGGCCAAGGGACAATGGTCACCGTCTCTTCAG"
    #         # (GGGGS)3 linker (nucleotide)
    #         "GGTGGCGGTTCAGGCGGAGGTGGCTCTGGCGGTGGCGGATCG"
    #         # VL domain (nucleotide)
    #         "GACATCCAGATGACCCAGTCTCCATCCTCCCTGTCTGCATCTGTAGGAGACAGAGTCACCATCACTTGCCGGGCAAGTCAGGGCATTAGCAGTTGGCTGGCCTGGTATCAGCAGAAACCAGGGAAAGCCCCTAAGCTCCTGATCTATGCTGCATCCAGTTTGCAAAGTGGGGTCCCATCAAGGTTCAGTGGCAGTGGATCTGGGACAGATTTCACTCTCACCATCAGCAGTCTGCAACCTGAAGATTTTGCAACTTACTACTGTCAGCAAGCTAACAGCTTCCCTTATACGTTCGGCCAAGGGACCAAGGTGGAAATCAAACGG"
    #     )

    nt_sample_result = nt_numbering.run_on_sequence("header", NT_QUERY, Scheme.IMGT)

    if isinstance(nt_sample_result, list):
        print(json.dumps([asdict(result) for result in nt_sample_result], indent=4))
    else:
        print(json.dumps(asdict(nt_sample_result), indent=4))

    # aa_numbering = create_riot_aa(return_all_domains=True)
    # AA_QUERY = "MPSSAVGVLGEAWYSLGGPDSSCAASGFTFSSYAMSWVRQAPGKGLEWVSSIANKGHETRYVDSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKYAGTFDYWGQGTLVTVSSGGGGSGGGGSGGGGSTDIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASMLQSGVPSRFSGSGSGTDFTLTISSLQPEYFATYYCQQARSWPPTFGQGDQGGNQTGRPHIIIAITGATHHHHHHGAAEQKLISEEDLNGAA"

    # aa_sample_result = aa_numbering.run_on_sequence("header", AA_QUERY, Scheme.IMGT)
    # # Handle both single result and list result
    # if isinstance(aa_sample_result, list):
    #     print(json.dumps([asdict(result) for result in aa_sample_result], indent=4))
    # else:
    #     print(json.dumps(asdict(aa_sample_result), indent=4))
