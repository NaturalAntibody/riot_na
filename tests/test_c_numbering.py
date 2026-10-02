import csv
from pathlib import Path

import pytest
from click.testing import CliRunner

from riot_na import ChainType, CScheme, Organism, Scheme
from riot_na.alignment.translate import translate
from riot_na.api.riot_numbering import (
    RiotNumberingAA,
    RiotNumberingNT,
    get_or_create_riot_aa,
    get_or_create_riot_nt,
)
from riot_na.cli import run_riot
from riot_na.common.fasta import parse_fasta
from riot_na.common.io import write_airr_iter_to_csv
from riot_na.config import GENE_DB_DIR
from riot_na.data.model import AirrRearrangementEntryAA
from riot_na.data.scheme_definitions import get_c_scheme_positions
from riot_na.data.scheme_mapping_facade import CSchemeMappingFacade
from tests.test_e2e import SEQUENCES_AA

IMGT_TABLE_PATH = Path(__file__).parents[1] / "notebooks" / "data_processing" / "imgt_c_numbering.tsv"

VH = (
    "QVQLLESGGGLVQPGGSLRLSCAASGFTFSTAAMSWVRQAPGKGLEWVSGISGSGSSTYYADSVKGRFTISRDNSKNTLYLQMNSLRATAVYYCARELSYLYSGYYF"
    "DYWGQGTLVTVSS"
)
VK = "EIVLTQSPATLSLSPGERATLSCRASKGVSTSGYSYLHWYQQKPGQAPRLLIYLASYLESGVPARFSGSGSGTDFTLTISSLEPEDFAVYYCQHSRDLPLTFGGGTKVEIK"
VL = "QSALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYDVSNRPSGVSNRFSGSKSGNTASLTISGLQAEDEADYYCSSYTSSSTLVVFGGGTKLTVL"
LIGHT_GENES = ("IGKC", "IGLC1", "IGLC2", "IGLC3", "IGLC6", "IGLC7")


def read_germlines(path: Path) -> dict[str, str]:
    return {record.description.split("\t")[0]: record.seq for record in parse_fasta(path)}


C_GERMLINES_AA = {
    gene: sequence
    for locus in ("igh", "igk", "igl")
    for gene, sequence in read_germlines(
        GENE_DB_DIR / "gene_db" / "aa_genes" / "c_genes" / "human" / f"{locus}.fasta"
    ).items()
}


def secreted_ighg(gene: str) -> str:
    # germlines in the database are membrane-bound forms, replace membrane exons with CHS
    germline = C_GERMLINES_AA[gene]
    return germline[: germline.index("ELQLEES")] + "GK"


@pytest.fixture(name="riot_aa")
def fixture_riot_aa() -> RiotNumberingAA:
    return get_or_create_riot_aa()


@pytest.fixture(name="riot_nt")
def fixture_riot_nt() -> RiotNumberingNT:
    return get_or_create_riot_nt()


def assert_positional_mapping_consistent(query: str, airr: AirrRearrangementEntryAA):
    assert airr.c_scheme_residue_mapping is not None
    assert airr.c_positional_scheme_mapping is not None
    assert list(airr.c_positional_scheme_mapping.values()) == list(airr.c_scheme_residue_mapping)
    for position, scheme_position in airr.c_positional_scheme_mapping.items():
        assert query[position] == airr.c_scheme_residue_mapping[scheme_position]


def imgt_table_numbering(gene: str, c_scheme: CScheme) -> dict[str, str]:
    with IMGT_TABLE_PATH.open() as table_file:
        rows = [row for row in csv.DictReader(table_file, delimiter="\t") if row["gene"] == gene]

    match c_scheme:
        case CScheme.EU:
            return {row["eu"]: row["residue"] for row in rows if row["eu"] != "-"}
        case CScheme.KABAT:
            # Kabat 107A is written as 107.1
            return {row["kabat"].replace("A", ".1"): row["residue"] for row in rows}
        case CScheme.IMGT:
            return {f"{row['region']}_{row['imgt_position']}": row["residue"] for row in rows}


@pytest.mark.parametrize("c_scheme", list(CScheme))
def test_ighg1_numbering_follows_imgt_table(riot_aa: RiotNumberingAA, c_scheme: CScheme):
    query = VH + secreted_ighg("IGHG1")
    airr = riot_aa.run_on_sequence("", query, Scheme.KABAT, c_scheme=c_scheme)

    assert airr.c_call == "IGHG1"
    assert airr.c_numbering_scheme == c_scheme.value
    # germline is membrane-bound, so CHS (EU 446-447) is not aligned
    assert airr.c_scheme_residue_mapping == dict(list(imgt_table_numbering("IGHG1", c_scheme).items())[:-2])
    assert_positional_mapping_consistent(query, airr)
    # V numbering is not affected
    assert airr.numbering_scheme == Scheme.KABAT.value
    assert airr.scheme_residue_mapping and "113" in airr.scheme_residue_mapping


@pytest.mark.parametrize("c_scheme", list(CScheme))
@pytest.mark.parametrize(
    "gene, query", [("IGKC", VK + C_GERMLINES_AA["IGKC"]), ("IGLC1", VL + C_GERMLINES_AA["IGLC1"])]
)
def test_light_chain_numbering_follows_imgt_table(riot_aa: RiotNumberingAA, gene: str, query: str, c_scheme: CScheme):
    airr = riot_aa.run_on_sequence("", query, Scheme.KABAT, c_scheme=c_scheme)
    expected = imgt_table_numbering(gene, c_scheme)

    assert airr.c_call == gene
    if not expected:
        # EU numbering is not defined for lambda chains
        assert airr.c_numbering_scheme is None
        assert airr.c_scheme_residue_mapping is None
        return

    assert airr.c_numbering_scheme == c_scheme.value
    assert airr.c_scheme_residue_mapping == expected
    assert_positional_mapping_consistent(query, airr)
    # Kabat numbering continues the variable domain one
    assert airr.scheme_residue_mapping and list(airr.scheme_residue_mapping)[-1] == "107"


@pytest.mark.parametrize("gene", LIGHT_GENES[2:])
def test_lambda_numbering(riot_aa: RiotNumberingAA, gene: str):
    # IGLC genes take IGLC1 positions
    query = VL + C_GERMLINES_AA[gene]
    for c_scheme in (CScheme.KABAT, CScheme.IMGT):
        airr = riot_aa.run_on_sequence("", query, c_scheme=c_scheme)
        assert airr.c_call == gene
        assert airr.c_scheme_residue_mapping
        assert list(airr.c_scheme_residue_mapping.values()) == list(C_GERMLINES_AA[gene])
        assert list(airr.c_scheme_residue_mapping) == list(imgt_table_numbering("IGLC1", c_scheme))


@pytest.mark.parametrize(
    "gene, c_scheme, expected, absent",
    [
        # IGHG2 lacks G236 (EU)
        (
            "IGHG2",
            CScheme.EU,
            {"219": "C", "220": "C", "222": "V", "233": "P", "234": "V", "235": "A", "237": "G"},
            ["221", "236"],
        ),
        ("IGHG2", CScheme.KABAT, {"246": "P", "247": "V", "248": "A", "250": "G"}, ["249"]),
        ("IGHG2", CScheme.IMGT, {"H_12": "P", "CH2_1.5": "P", "CH2_1.3": "P", "CH2_1.1": "A"}, ["CH2_1.4", "H_13"]),
        # IGHG3 hinge has 47 residues without EU positions
        (
            "IGHG3",
            CScheme.EU,
            {"218": "K", "218.1": "T", "218.2": "P", "219": "L", "228": "R", "228.1": "C", "228.45": "R", "229": "C"},
            ["228.46"],
        ),
        (
            "IGHG3",
            CScheme.KABAT,
            {"229": "T", "230": "P", "241": "R", "241.1": "C", "241.45": "R", "242": "C"},
            ["231"],
        ),
        ("IGHG3", CScheme.IMGT, {"H1_1": "E", "H1_17": "P", "H2_1": "E", "H4_15": "P", "CH2_1.6": "A"}, ["H_1"]),
        (
            "IGHG4",
            CScheme.EU,
            {"220": "G", "224": "P", "228": "S", "234": "F", "235": "L", "409": "R", "445": "L"},
            ["221", "222", "223"],
        ),
        ("IGHG4", CScheme.KABAT, {"230": "G", "237": "P", "241": "S", "247": "F"}, ["231", "232", "236"]),
        ("IGHG4", CScheme.IMGT, {"H_10": "S", "H_12": "P", "CH2_1.3": "F", "CH2_84.4": "N"}, ["H_13"]),
    ],
)
def test_ighg_numbering(
    riot_aa: RiotNumberingAA, gene: str, c_scheme: CScheme, expected: dict[str, str], absent: list[str]
):
    query = VH + secreted_ighg(gene)
    airr = riot_aa.run_on_sequence("", query, Scheme.IMGT, c_scheme=c_scheme)

    assert airr.c_call == gene
    assert airr.c_scheme_residue_mapping is not None
    for scheme_position, residue in expected.items():
        assert airr.c_scheme_residue_mapping[scheme_position] == residue, scheme_position
    for scheme_position in absent:
        assert scheme_position not in airr.c_scheme_residue_mapping
    assert_positional_mapping_consistent(query, airr)


def test_c_numbering_insertion_and_deletion(riot_aa: RiotNumberingAA):
    constant = secreted_ighg("IGHG1")
    n297 = constant.index("QYNSTY") + 2

    with_insertion = VH + constant[: n297 + 1] + "W" + constant[n297 + 1 :]
    eu = riot_aa.run_on_sequence("", with_insertion, c_scheme=CScheme.EU).c_scheme_residue_mapping
    assert eu and (eu["297"], eu["297.1"], eu["298"]) == ("N", "W", "S")
    # insertions take unused scheme positions next to them if there are any
    kabat = riot_aa.run_on_sequence("", with_insertion, c_scheme=CScheme.KABAT).c_scheme_residue_mapping
    assert kabat and (kabat["314"], kabat["315"], kabat["317"]) == ("N", "W", "S")
    imgt = riot_aa.run_on_sequence("", with_insertion, c_scheme=CScheme.IMGT).c_scheme_residue_mapping
    assert imgt and (imgt["CH2_84.4"], imgt["CH2_84.5"], imgt["CH2_85.4"]) == ("N", "W", "S")

    with_deletion = VH + constant[:n297] + constant[n297 + 1 :]
    airr = riot_aa.run_on_sequence("", with_deletion, c_scheme=CScheme.EU)
    assert airr.c_scheme_residue_mapping
    assert "297" not in airr.c_scheme_residue_mapping
    assert (airr.c_scheme_residue_mapping["296"], airr.c_scheme_residue_mapping["298"]) == ("Y", "S")
    assert_positional_mapping_consistent(with_deletion, airr)


def test_c_numbering_partial_constant_region(riot_aa: RiotNumberingAA):
    # CH1 without first 3 residues
    query = VH + secreted_ighg("IGHG1")[3:120]
    airr = riot_aa.run_on_sequence("", query, c_scheme=CScheme.EU)
    assert airr.c_scheme_residue_mapping
    assert list(airr.c_scheme_residue_mapping)[0] == "121"
    assert_positional_mapping_consistent(query, airr)


def test_membrane_exons_are_not_numbered(riot_aa: RiotNumberingAA):
    query = VH + C_GERMLINES_AA["IGHG1"]
    airr = riot_aa.run_on_sequence("", query, c_scheme=CScheme.EU)
    assert airr.c_scheme_residue_mapping
    assert list(airr.c_scheme_residue_mapping)[-1] == "445"
    assert_positional_mapping_consistent(query, airr)


def test_no_c_numbering(riot_aa: RiotNumberingAA):
    airr = riot_aa.run_on_sequence("", VH + secreted_ighg("IGHG1"))
    assert airr.c_call == "IGHG1"
    assert airr.c_numbering_scheme is None
    assert airr.c_scheme_residue_mapping is None
    assert airr.c_positional_scheme_mapping is None

    # only IGHG1-4, IGKC and IGLC have C scheme mappings
    airr = riot_aa.run_on_sequence("", VH + C_GERMLINES_AA["IGHM"], c_scheme=CScheme.EU)
    assert airr.c_call == "IGHM"
    assert airr.c_scheme_residue_mapping is None


def test_c_numbering_multiple_domains():
    riot_aa = get_or_create_riot_aa(return_all_domains=True)
    query = VH + secreted_ighg("IGHG1") + "GGGGSGGGGS" + SEQUENCES_AA["L"]
    heavy, light = riot_aa.run_on_sequence("", query, return_all_domains=True, c_scheme=CScheme.EU)

    assert heavy.c_scheme_residue_mapping and heavy.c_scheme_residue_mapping["297"] == "N"
    assert light.c_scheme_residue_mapping and light.c_scheme_residue_mapping["214"] == "C"


@pytest.mark.parametrize("c_scheme", list(CScheme))
@pytest.mark.parametrize(
    "v_gene, cdr3, j_gene, j_frame, c_gene",
    [
        ("IGHV3-23*01", "GATCGGGGGTAC", "IGHJ4*02", 2, "IGHG1"),
        ("IGKV1-27*01", "CAGCAGAGTTACAGTACCCCT", "IGKJ1*01", 1, "IGKC"),
        ("IGLV2-14*01", "AGCTCATATACAAGCAGCAGCACT", "IGLJ2*01", 1, "IGLC2"),
    ],
)
def test_c_numbering_nucleotides(
    riot_nt: RiotNumberingNT,
    riot_aa: RiotNumberingAA,
    v_gene: str,
    cdr3: str,
    j_gene: str,
    j_frame: int,
    c_gene: str,
    c_scheme: CScheme,
):
    nt_germlines = read_germlines(GENE_DB_DIR / "gene_db" / "v_genes" / "human.fasta")
    for locus in ("igh", "igk", "igl"):
        nt_germlines |= read_germlines(GENE_DB_DIR / "gene_db" / "j_genes" / "human" / f"{locus}.fasta")
        nt_germlines |= read_germlines(GENE_DB_DIR / "gene_db" / "c_genes" / "human" / f"{locus}.fasta")

    v_sequence = nt_germlines[v_gene]
    c_sequence = nt_germlines[c_gene]
    # in frame V, CDR3, J without the splice nucleotide and C without stop codon
    query = (
        v_sequence[: len(v_sequence) // 3 * 3]
        + cdr3
        + nt_germlines[j_gene][j_frame:-1]
        + (c_sequence[:-3] if translate(c_sequence[-3:], 0) == "*" else c_sequence)
    )

    airr_nt = riot_nt.run_on_sequence("", query, c_scheme=c_scheme)
    airr_aa = riot_aa.run_on_sequence("", translate(query, 0), c_scheme=c_scheme)

    assert airr_nt.c_call == c_gene
    assert airr_nt.c_scheme_residue_mapping == airr_aa.c_scheme_residue_mapping
    assert airr_nt.c_positional_scheme_mapping == airr_aa.c_positional_scheme_mapping
    if c_scheme == CScheme.EU and c_gene.startswith("IGLC"):
        assert airr_nt.c_scheme_residue_mapping is None
    else:
        assert airr_nt.c_numbering_scheme == c_scheme.value
        assert airr_nt.c_scheme_residue_mapping


@pytest.mark.parametrize("c_scheme", list(CScheme))
def test_c_scheme_mappings(c_scheme: CScheme):
    facade = CSchemeMappingFacade(c_scheme)

    for gene in ("IGHG1", "IGHG2", "IGHG3", "IGHG4", *LIGHT_GENES):
        mapping = facade.get_mapping(Organism.HOMO_SAPIENS, gene)
        if c_scheme == CScheme.EU and gene.startswith("IGLC"):
            assert mapping is None
            continue

        chain_type = ChainType.HEAVY if gene.startswith("IGH") else ChainType.LIGHT
        germline = secreted_ighg(gene)[:-2] if chain_type == ChainType.HEAVY else C_GERMLINES_AA[gene]
        assert mapping
        assert set(mapping) <= {"M", "D", "I"}
        assert mapping.count("M") + mapping.count("D") <= len(get_c_scheme_positions(c_scheme, chain_type))
        assert mapping.count("M") + mapping.count("I") == len(germline)
        # mouse sequences are aligned to human C germlines
        assert facade.get_mapping(Organism.MUS_MUSCULUS, gene) == mapping

    assert facade.get_mapping(Organism.HOMO_SAPIENS, "IGHM") is None


def test_c_numbering_csv_output(riot_aa: RiotNumberingAA, tmp_path: Path):
    airr = riot_aa.run_on_sequence("", VH + secreted_ighg("IGHG1"), c_scheme=CScheme.IMGT)
    output_path = tmp_path / "output.csv"
    write_airr_iter_to_csv(output_path, AirrRearrangementEntryAA, [airr])

    with output_path.open() as output_file:
        row = next(csv.DictReader(output_file))
    assert row["c_numbering_scheme"] == "imgt"
    assert row["c_scheme_residue_mapping"]
    assert row["c_positional_scheme_mapping"]


@pytest.mark.parametrize("multiple_domains", [False, True])
def test_cli_single_sequence_csv_output(tmp_path: Path, multiple_domains: bool):
    output_path = tmp_path / "output.csv"
    args = ["-s", VH + secreted_ighg("IGHG1"), "--input-type", "aa", "--species", "human", "--c-scheme", "eu"]
    args += ["-o", str(output_path), "--multiple-domains", str(multiple_domains)]

    result = CliRunner().invoke(run_riot, args)

    assert result.exit_code == 0, result.output
    with output_path.open() as output_file:
        rows = list(csv.DictReader(output_file))
    assert len(rows) == 1
    assert rows[0]["c_call"] == "IGHG1"
    assert rows[0]["c_numbering_scheme"] == "eu"
