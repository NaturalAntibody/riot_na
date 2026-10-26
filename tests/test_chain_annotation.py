import re

import pytest

from riot_na import ChainType, CScheme, Organism
from riot_na.alignment.translate import translate
from riot_na.api.chain_annotation import reverse_complement
from riot_na.api.riot_numbering import (
    RiotNumberingAA,
    RiotNumberingNT,
    get_or_create_riot_aa,
    get_or_create_riot_nt,
    v_domain_span_aa,
    v_domain_span_nt,
)
from riot_na.data.model import ConstantDomain
from riot_na.data.scheme_mapping_facade import CSchemeMappingFacade
from riot_na.schemes.constant_domains import GermlineDomain, get_germline_domains
from tests.test_c_numbering import (
    C_GERMLINES_AA,
    NT_LINKER,
    VH,
    VK,
    VL,
    nt_rearrangement,
    secreted_ighg,
)

IGHG1 = secreted_ighg("IGHG1")
FC = IGHG1[IGHG1.index("DKTHT") :]
LINKER = "GGGGSGGGGSGGGGS"
VHH = (
    "QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAGRYGSSWY"
    "PDSYDYWGQGTQVTVSS"
)
# CAR CD8a hinge and transmembrane domain, CD28 signalling domain
CAR_TAIL = (
    "TTTPAPRPPTPAPTIASQPLSLRPEACRPAAGGAVHTRGLDFACDIYIWAPLAGTCGVLLLSLVITLYCKRGRKKLLYIFKQPFMRPVQTTQEEDGCSCRFPEEEEGGCEL"
)
VK_NT = nt_rearrangement("IGKV1-27*01", "CAGCAGAGTTACAGTACCCCT", "IGKJ1*01", 1)


@pytest.fixture(name="riot_aa")
def fixture_riot_aa() -> RiotNumberingAA:
    return get_or_create_riot_aa(return_all_domains=True)


@pytest.fixture(name="riot_nt")
def fixture_riot_nt() -> RiotNumberingNT:
    return get_or_create_riot_nt(return_all_domains=True)


def eu_range(constant_domain: ConstantDomain) -> tuple[str, str]:
    assert constant_domain.c_scheme_residue_mapping
    positions = list(constant_domain.c_scheme_residue_mapping)
    return positions[0], positions[-1]


def test_germline_domains():
    imgt_mapping = CSchemeMappingFacade(CScheme.IMGT).get_mapping(Organism.HOMO_SAPIENS, "IGHG1")
    assert imgt_mapping

    assert get_germline_domains(imgt_mapping, ChainType.HEAVY, len(C_GERMLINES_AA["IGHG1"])) == [
        GermlineDomain("CH1", 0, 98),
        GermlineDomain("H", 98, 113),
        GermlineDomain("CH2", 113, 223),
        GermlineDomain("CH3", 223, 328),
        # membrane exons
        GermlineDomain("M", 328, len(C_GERMLINES_AA["IGHG1"])),
    ]


def test_heavy_chain(riot_aa: RiotNumberingAA):
    query = VH + IGHG1
    chain = riot_aa.annotate_chain("heavy", query, c_scheme=CScheme.EU)

    assert chain.layout == "VH-CH1-H-CH2-CH3"
    assert [rearrangement.locus for rearrangement in chain.rearrangements] == ["igh"]
    assert [
        (domain.domain, domain.c_call, domain.rearrangement_index, eu_range(domain))
        for domain in chain.constant_domains
    ] == [
        ("CH1", "IGHG1", 0, ("118", "215")),
        ("H", "IGHG1", 0, ("216", "230")),
        ("CH2", "IGHG1", 0, ("231", "340")),
        ("CH3", "IGHG1", 0, ("341", "445")),
    ]
    # domains cover the constant region, CH3 includes CHS (GK) of secreted IgG, which is not numbered
    assert "".join(domain.sequence_aa for domain in chain.constant_domains) == IGHG1
    assert (chain.constant_domains[0].sequence_start, chain.constant_domains[-1].sequence_end) == (
        len(VH) + 1,
        len(query),
    )
    for domain in chain.constant_domains:
        assert domain.c_positional_scheme_mapping
        for position, scheme_position in domain.c_positional_scheme_mapping.items():
            assert query[position] == domain.c_scheme_residue_mapping[scheme_position]  # type: ignore[index]


def test_rearrangements_are_multiple_domains_results(riot_aa: RiotNumberingAA):
    query = VH + IGHG1 + LINKER + VK
    chain = riot_aa.annotate_chain("", query, c_scheme=CScheme.EU)
    domains = riot_aa.run_on_sequence("", query, return_all_domains=True, c_scheme=CScheme.EU)

    assert chain.rearrangements == domains
    assert chain.layout == "VH-CH1-H-CH2-CH3-[15]-VK"


@pytest.mark.parametrize(
    "query, layout, domains",
    [
        # CrossMab chains
        (VK + IGHG1, "VK-CH1-H-CH2-CH3", [("CH1", "IGHG1"), ("H", "IGHG1"), ("CH2", "IGHG1"), ("CH3", "IGHG1")]),
        (VH + C_GERMLINES_AA["IGKC"], "VH-CL", [("CL", "IGKC")]),
        # scFv-Fc
        (VH + LINKER + VK + FC, "VH-[15]-VK-H-CH2-CH3", [("H", "IGHG1"), ("CH2", "IGHG1"), ("CH3", "IGHG1")]),
        # Fc chain of a knobs-into-holes bispecific
        (FC, "H-CH2-CH3", [("H", "IGHG1"), ("CH2", "IGHG1"), ("CH3", "IGHG1")]),
        # VHH with human Fc
        (VHH + FC, "VH-H-CH2-CH3", [("H", "IGHG1"), ("CH2", "IGHG1"), ("CH3", "IGHG1")]),
        # Fab heavy chain with CH1 fragment
        (VH + "ASTKGPS", "VH-CH1", [("CH1", "IGHG1")]),
        # constant domains of other loci and genes without domain definitions
        (VL + C_GERMLINES_AA["IGLC2"], "VL-CL", [("CL", "IGLC2")]),
        (VH + C_GERMLINES_AA["IGHM"], "VH-C", [("C", "IGHM")]),
    ],
    ids=["VK-CH1", "VH-CL", "scFv-Fc", "Fc", "VHH-Fc", "Fab", "lambda", "IgM"],
)
def test_constant_domains(riot_aa: RiotNumberingAA, query: str, layout: str, domains: list[tuple[str, str]]):
    chain = riot_aa.annotate_chain("", query)

    assert chain.layout == layout
    assert [(domain.domain, domain.c_call) for domain in chain.constant_domains] == domains


def test_hybrid_constant_region(riot_aa: RiotNumberingAA):
    # eculizumab-like IgG2 CH1 and hinge with IgG4 CH2 and CH3
    hybrid = secreted_ighg("IGHG2")[:110] + secreted_ighg("IGHG4")[110:]
    chain = riot_aa.annotate_chain("", VH + hybrid, c_scheme=CScheme.EU)

    assert [(domain.domain, domain.c_call) for domain in chain.constant_domains] == [
        ("CH1", "IGHG2"),
        ("H", "IGHG2"),
        ("CH2", "IGHG4"),
        ("CH3", "IGHG4"),
    ]
    # hinge is numbered as IGHG2 hinge, CH2 as IGHG4 CH2
    hinge, ch2 = chain.constant_domains[1].c_scheme_residue_mapping, chain.constant_domains[2].c_scheme_residue_mapping
    assert hinge and ch2
    assert (hinge["219"], hinge["220"], "221" in hinge) == ("C", "C", False)
    assert (ch2["234"], ch2["235"]) == ("F", "L")


def test_repeated_constant_regions(riot_aa: RiotNumberingAA):
    # single chain Fc following a scFv, both Fc regions are numbered
    query = VH + LINKER + VK + "GGGG" + FC + LINKER + "GGGGS" + FC
    chain = riot_aa.annotate_chain("", query, c_scheme=CScheme.EU)

    assert chain.layout == "VH-[15]-VK-[4]-H-CH2-CH3-[20]-H-CH2-CH3"
    assert [(domain.domain, eu_range(domain)) for domain in chain.constant_domains] == 2 * [
        ("H", ("221", "230")),
        ("CH2", ("231", "340")),
        ("CH3", ("341", "445")),
    ]
    # constant domains follow the VK domain
    assert {domain.rearrangement_index for domain in chain.constant_domains} == {1}


def test_constant_domains_preceding_v_domains(riot_aa: RiotNumberingAA):
    chain = riot_aa.annotate_chain("", FC + LINKER + VHH)

    assert chain.layout == "H-CH2-CH3-[15]-VH"
    assert {domain.rearrangement_index for domain in chain.constant_domains} == {None}


def test_constant_region_fragments_follow_v_domains(riot_aa: RiotNumberingAA):
    ck_start = C_GERMLINES_AA["IGKC"][:12]

    # CL fragment of the J-C junction (CL derived linker), and in front of the V domain
    assert riot_aa.annotate_chain("", VK + ck_start).layout == "VK-CL"
    assert riot_aa.annotate_chain("", ck_start + LINKER + VK).layout == "[27]-VK"
    assert riot_aa.annotate_chain("", VK + LINKER + ck_start).layout == "VK-[27]"

    # CL derived linker of DVD-Ig light chains, without the J-C junction residue (R108)
    assert riot_aa.annotate_chain("", VK + "TVAAP" + VL).layout == "VK-CL-VL"


def test_constant_domains_are_aligned_to_c_genes_of_v_domain_species_first(riot_aa: RiotNumberingAA):
    # C genes of the species of the preceding V domain win ties, in the same order in every process
    assert riot_aa.constant_region_aligner
    assert riot_aa.constant_region_aligner.aligned_species(Organism.VICUGNA_PACOS) == [
        Organism.VICUGNA_PACOS,
        Organism.HOMO_SAPIENS,
    ]
    assert riot_aa.constant_region_aligner.aligned_species(Organism.HOMO_SAPIENS) == [Organism.HOMO_SAPIENS]


def test_tandem_fabs(riot_aa: RiotNumberingAA):
    query = VH + IGHG1[:98] + LINKER + VH + IGHG1
    chain = riot_aa.annotate_chain("", query)

    assert chain.layout == "VH-CH1-[15]-VH-CH1-H-CH2-CH3"
    assert [domain.rearrangement_index for domain in chain.constant_domains] == [0, 1, 1, 1, 1]


@pytest.mark.parametrize("query", [VH + LINKER + VK + CAR_TAIL, VK + LINKER + VH + CAR_TAIL], ids=["VH-VK", "VK-VH"])
def test_non_antibody_regions_are_not_constant_domains(riot_aa: RiotNumberingAA, query: str):
    chain = riot_aa.annotate_chain("", query)
    assert not chain.constant_domains


def test_nucleotide_chain(riot_nt: RiotNumberingNT, riot_aa: RiotNumberingAA):
    heavy_chain = nt_rearrangement("IGHV3-23*01", "GGCTAC", "IGHJ4*02", 2, "IGHG1")

    chain = riot_nt.annotate_chain("", heavy_chain, c_scheme=CScheme.EU)
    assert chain.layout == "VH-CH1-H-CH2-CH3"
    assert not chain.rev_comp
    assert chain.constant_domains[-1].sequence_end == len(heavy_chain)

    # numbering is the same as of the translated chain, positions refer to the first codon nucleotides
    aa_chain = riot_aa.annotate_chain(
        "", "".join(domain.sequence_aa for domain in chain.constant_domains), c_scheme=CScheme.EU
    )
    assert [domain.c_scheme_residue_mapping for domain in chain.constant_domains] == [
        domain.c_scheme_residue_mapping for domain in aa_chain.constant_domains
    ]
    first_domain = chain.constant_domains[0]
    assert first_domain.c_positional_scheme_mapping
    assert list(first_domain.c_positional_scheme_mapping)[:2] == [
        first_domain.sequence_start - 1,
        first_domain.sequence_start + 2,
    ]

    # reverse strand
    reverse_chain = riot_nt.annotate_chain("", reverse_complement(heavy_chain), c_scheme=CScheme.EU)
    assert reverse_chain.rev_comp
    assert reverse_chain.sequence == heavy_chain
    assert reverse_chain.constant_domains == chain.constant_domains


def test_nucleotide_chain_without_v_domains(riot_nt: RiotNumberingNT):
    heavy_chain = nt_rearrangement("IGHV3-23*01", "GGCTAC", "IGHJ4*02", 2, "IGHG1")
    fc = heavy_chain[heavy_chain.index("GACAAAACTCACACATGC") :]

    for query, rev_comp in ((fc, False), (reverse_complement(fc), True)):
        chain = riot_nt.annotate_chain("", query)
        assert chain.rev_comp == rev_comp
        assert chain.layout == "H-CH2-CH3"


def test_nucleotide_chain_with_multiple_domains(riot_nt: RiotNumberingNT):
    # IgG with scFv fused to the C-terminus of the heavy chain
    query = (
        nt_rearrangement("IGHV3-23*01", "GGCTAC", "IGHJ4*02", 2, "IGHG1")
        + NT_LINKER
        + nt_rearrangement("IGHV1-69*01", "GATCGGGGCTAC", "IGHJ4*02", 2)
        + NT_LINKER
        + nt_rearrangement("IGKV3-11*01", "TGGCCGCTC", "IGKJ4*01", 1)
    )
    chain = riot_nt.annotate_chain("", query)

    # unannotated nucleotides of linkers
    assert re.fullmatch(r"VH-CH1-H-CH2-CH3-\[4\d\]-VH-\[4\d\]-VK", chain.layout)


def test_v_domains_without_j_end_with_fr4(riot_aa: RiotNumberingAA, riot_nt: RiotNumberingNT):
    # V domains end with FR4 following the J motif (W-G-x-G), if J is not found
    aa_alignments = riot_aa.vjc_aligner_aa.produce_aa_alignments(VH + "GGGGS" + VK)[0]
    assert v_domain_span_aa(aa_alignments).end == len(VH)
    aa_alignments.j = None
    assert v_domain_span_aa(aa_alignments).end == len(VH)

    # J germline ends with the nucleotide of the first C codon, which matches the first linker nucleotide (G)
    heavy_chain = nt_rearrangement("IGHV3-23*01", "GGCTAC", "IGHJ4*02", 2)
    nt_alignments = riot_nt.vdjc_aligner.produce_nt_alignments(heavy_chain + "GGTGGCGGTGGCTCG" + VK_NT)[0]
    assert v_domain_span_nt(nt_alignments).end == len(heavy_chain)
    nt_alignments.j = None
    assert v_domain_span_nt(nt_alignments).end == len(heavy_chain)


def test_nucleotide_v_domains_end_with_fr4_residues(riot_nt: RiotNumberingNT):
    # synonymous codons at the FR4 end (TCC TCA -> AGC AGT), nucleotide J alignment ends before them
    heavy_chain = nt_rearrangement("IGHV3-23*01", "GGCTAC", "IGHJ4*02", 2, "IGHG1")
    fr4_end = heavy_chain.index("GCCTCCACCAAGGGC")
    recoded_chain = heavy_chain[: fr4_end - 6] + "AGCAGT" + heavy_chain[fr4_end:]
    assert translate(recoded_chain, 0) == translate(heavy_chain, 0)

    assert riot_nt.annotate_chain("", recoded_chain).layout == "VH-CH1-H-CH2-CH3"


def test_chain_annotation_requires_multiple_domains_numbering():
    with pytest.raises(ValueError, match="annotate_chain"):
        get_or_create_riot_aa().annotate_chain("", VH + IGHG1)
    with pytest.raises(ValueError, match="annotate_chain"):
        get_or_create_riot_nt().annotate_chain("", "GAGGTGCAGCTGGTGGAGTCTGGG")
