import csv
import json
from pathlib import Path
from typing import Optional

import pytest

from riot_na import CScheme, MoleculeAnnotation
from riot_na.api.chain_annotation import reverse_complement
from riot_na.api.molecule_annotation import (
    _FabArm,
    _one_to_one_partners,
    _VDomain,
    binding_site_rows,
    molecule_layout_rows,
)
from riot_na.api.riot_numbering import (
    RiotNumberingAA,
    get_or_create_riot_aa,
    get_or_create_riot_nt,
)
from riot_na.common.io import (
    read_v_domain_pairs,
    write_binding_sites_csv,
    write_molecule_layout_csv,
    write_molecules_json,
)
from riot_na.data.model import ChainDomain
from tests.test_c_numbering import (
    C_GERMLINES_AA,
    VH,
    VK,
    VL,
    nt_rearrangement,
    secreted_ighg,
)

VH2 = (
    "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYWMHWVRQAPGQGLEWMGEINPSNGRTNYNEKFKSRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARRDYDYDGRYF"
    "DYWGQGTTVTVSS"
)
VK2 = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
VHH = (
    "QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAGRYGSSWY"
    "PDSYDYWGQGTQVTVSS"
)
IGHG1 = secreted_ighg("IGHG1")
CH1 = IGHG1[:98]
FC = IGHG1[IGHG1.index("DKTHT") :]
CK = C_GERMLINES_AA["IGKC"]
LINKER = "GGGGSGGGGSGGGGS"


@pytest.fixture(name="riot_aa")
def fixture_riot_aa() -> RiotNumberingAA:
    return get_or_create_riot_aa(return_all_domains=True)


def binding_sites(molecule: MoleculeAnnotation) -> list[tuple[str, Optional[str], Optional[str], Optional[str]]]:
    """(type, heavy V domain, light V domain, candidates) of binding sites, V domains as chain:domain index"""
    return [
        (
            row["binding_site_type"],
            row["heavy_chain_id"] and f"{row['heavy_chain_id']}:{row['heavy_domain_index']}",
            row["light_chain_id"] and f"{row['light_chain_id']}:{row['light_domain_index']}",
            row["candidates"],
        )
        for row in binding_site_rows(molecule)
    ]


@pytest.mark.parametrize(
    "chains, expected",
    [
        ({"HC": VH + IGHG1, "LC": VK + CK}, [("Fab", "HC:1", "LC:1", None)]),
        # common light chain
        (
            {"HC1": VH + IGHG1, "HC2": VH2 + IGHG1, "LC": VK + CK},
            [("Fab", "HC1:1", "LC:1", None), ("Fab", "HC2:1", "LC:1", None)],
        ),
        # CrossMab, VL-CH1 pairs with VH-CL
        (
            {"HC1": VH + IGHG1, "LC1": VK + CK, "HC2": VK2 + IGHG1, "LC2": VH2 + CK},
            [("Fab", "HC1:1", "LC1:1", None), ("Fab", "LC2:1", "HC2:1", None)],
        ),
        # DVD-Ig, outer V domains joined by CH1 and CL derived linkers (annotated as CH1 and CL fragments)
        (
            {"HC": VH + "ASTKGP" + VH2 + IGHG1, "LC": VK + "TVAAP" + VK2 + CK},
            [("Fv", "HC:1", "LC:1", None), ("Fab", "HC:3", "LC:3", None)],
        ),
        # FIT-Ig
        (
            {"LC-HC": VK + CK + VH2 + IGHG1, "Fd": VH + CH1, "LC": VK2 + CK},
            [("Fab", "Fd:1", "LC-HC:1", None), ("Fab", "LC-HC:3", "LC:1", None)],
        ),
        # IgG with scFv fused to the heavy chain
        (
            {"HC": VH + IGHG1 + LINKER + VH2 + LINKER + VK2, "LC": VK + CK},
            [("Fab", "HC:1", "LC:1", None), ("scFv", "HC:6", "HC:7", None)],
        ),
        # tandem scFv (BiTE)
        (
            {"BiTE": VK + LINKER + VH + "GGGGS" + VH2 + LINKER + VK2},
            [("scFv", "BiTE:2", "BiTE:1", None), ("scFv", "BiTE:3", "BiTE:4", None)],
        ),
        # diabodies: homodimer and heterodimer
        ({"chain": VH + "GGGGS" + VK}, [("diabody", "chain:1", "chain:2", None)]),
        (
            {"c1": VH + "GGGGS" + VK2, "c2": VH2 + "GGGGS" + VK},
            [("diabody", "c1:1", "c2:2", None), ("diabody", "c2:1", "c1:2", None)],
        ),
        # single chain Fab
        ({"scFab": VK + CK + 6 * "GGGGS" + VH + CH1}, [("Fab", "scFab:3", "scFab:1", None)]),
        # single chain Fabs pair within their chains
        (
            {"scFab1": VK + CK + 8 * "GGGGS" + VH + IGHG1, "scFab2": VK2 + CK + 8 * "GGGGS" + VH2 + IGHG1},
            [("Fab", "scFab1:3", "scFab1:1", None), ("Fab", "scFab2:3", "scFab2:1", None)],
        ),
        # constant regions of C genes without domain definitions
        ({"HC": VH + C_GERMLINES_AA["IGHM"], "LC": VK + CK}, [("Fab", "HC:1", "LC:1", None)]),
        ({"VHH-Fc": VHH + FC}, [("single domain", "VHH-Fc:1", None, None)]),
        # one-armed antibody with a Fc chain
        ({"HC": VH + IGHG1, "Fc": FC, "LC": VK + CK}, [("Fab", "HC:1", "LC:1", None)]),
        ({"HC": VH + IGHG1}, [("unpaired", "HC:1", None, None)]),
        # heavy and light chains of two antibodies can't be paired from sequences
        (
            {"HC1": VH + IGHG1, "HC2": VH2 + IGHG1, "LC1": VK + CK, "LC2": VK2 + CK},
            [("Fab", "HC1:1", None, "LC1:1;LC2:1"), ("Fab", "HC2:1", None, "LC1:1;LC2:1")],
        ),
    ],
    ids=[
        "IgG",
        "common-LC",
        "CrossMab",
        "DVD-Ig",
        "FIT-Ig",
        "IgG-scFv",
        "BiTE",
        "diabody",
        "heterodimeric-diabody",
        "scFab",
        "scFab-chains",
        "IgM",
        "VHH-Fc",
        "one-armed",
        "HC-only",
        "bispecific",
    ],
)
def test_binding_sites(riot_aa: RiotNumberingAA, chains: dict[str, str], expected: list):
    molecule = riot_aa.annotate_molecule("molecule", chains)
    assert binding_sites(molecule) == expected


def test_binding_site_details(riot_aa: RiotNumberingAA):
    molecule = riot_aa.annotate_molecule(
        "molecule", {"HC": VH + IGHG1 + LINKER + VH2 + LINKER + VK2, "LC": VK + CK}, c_scheme=CScheme.EU
    )

    fab, scfv = molecule.binding_sites
    assert (fab.ambiguous, fab.linker_length) == (False, None)
    assert scfv.linker_length == len(LINKER)
    # chains are annotated with headers molecule_id|chain_id, as AIRR records
    assert molecule.chains["HC"].sequence_header == "molecule|HC"
    assert {rearrangement.sequence_header for rearrangement in molecule.chains["HC"].rearrangements} == {"molecule|HC"}
    assert molecule.chains["HC"].constant_domains[0].c_scheme_residue_mapping


BISPECIFIC = {"HC1": VH + IGHG1, "HC2": VH2 + IGHG1, "LC1": VK + CK, "LC2": VK2 + CK}
TANDEM_FAB = {"HC": VH + CH1 + "GGGGS" + VH2 + IGHG1, "LC1": VK + CK, "LC2": VK2 + CK}
DIABODIES = {"A": VH + "GGGGS" + VK2, "B": VH2 + "GGGGS" + VK, "C": VH + "GGGGS" + VK}
DVD_IG = {"HC": VH + "ASTKGP" + VH2 + IGHG1, "LC": VK + "TVAAP" + VK2 + CK}


@pytest.mark.parametrize(
    "chains, expected",
    [
        (TANDEM_FAB, [("Fab", "HC:1", None, "LC1:1;LC2:1"), ("Fab", "HC:3", None, "LC1:1;LC2:1")]),
        ({"H": VH, "L": VK}, [("single domain", "H:1", None, None), ("unpaired", None, "L:1", None)]),
        (
            DIABODIES,
            [
                ("diabody", "A:1", "A:2", "B:1;B:2;C:1;C:2"),
                ("diabody", "B:1", "B:2", "A:1;A:2;C:1;C:2"),
                ("diabody", "C:1", "C:2", "A:1;A:2;B:1;B:2"),
            ],
        ),
    ],
    ids=["tandem-Fab", "V-domains", "diabodies"],
)
def test_ambiguous_binding_sites(riot_aa: RiotNumberingAA, chains: dict[str, str], expected: list):
    molecule = riot_aa.annotate_molecule("molecule", chains)
    assert binding_sites(molecule) == expected


def numbered_binding_sites(
    molecule: MoleculeAnnotation,
) -> list[tuple[str, Optional[str], Optional[str], Optional[int]]]:
    """(type, heavy V domain, light V domain, linker length) of binding sites, V domains as chain:V domain number"""
    return [
        (
            row["binding_site_type"],
            row["heavy_chain_id"] and f"{row['heavy_chain_id']}:{row['heavy_v_domain_number']}",
            row["light_chain_id"] and f"{row['light_chain_id']}:{row['light_v_domain_number']}",
            row["linker_length"],
        )
        for row in binding_site_rows(molecule)
    ]


@pytest.mark.parametrize(
    "chains, pairs, expected",
    [
        ({"HC": VH + IGHG1, "LC": VK + CK}, [(1, 1, 2, 1)], [("Fab", "HC:1", "LC:1", None)]),
        # pairs which can't be told from the sequences
        (BISPECIFIC, [(1, 1, 4, 1), (2, 1, 3, 1)], [("Fab", "HC1:1", "LC2:1", None), ("Fab", "HC2:1", "LC1:1", None)]),
        # chains given by chain ID
        (
            BISPECIFIC,
            [("HC1", 1, "LC2", 1), ("HC2", 1, "LC1", 1)],
            [("Fab", "HC1:1", "LC2:1", None), ("Fab", "HC2:1", "LC1:1", None)],
        ),
        (
            {"HC1": VH + IGHG1, "HC2": VH2 + IGHG1, "LC": VK + CK},
            [(1, 1, 3, 1), (2, 1, 3, 1)],
            [("Fab", "HC1:1", "LC:1", None), ("Fab", "HC2:1", "LC:1", None)],
        ),
        # common heavy chain with a kappa and a lambda light chain
        (
            {"HC": VH + IGHG1, "LC1": VK + CK, "LC2": VL + C_GERMLINES_AA["IGLC2"]},
            [(1, 1, 2, 1), (1, 1, 3, 1)],
            [("Fab", "HC:1", "LC1:1", None), ("Fab", "HC:1", "LC2:1", None)],
        ),
        (
            {"HC1": VH + IGHG1, "LC1": VK + CK, "HC2": VK2 + IGHG1, "LC2": VH2 + CK},
            [(1, 1, 2, 1), (3, 1, 4, 1)],
            [("Fab", "HC1:1", "LC1:1", None), ("Fab", "LC2:1", "HC2:1", None)],
        ),
        (TANDEM_FAB, [(1, 1, 2, 1), (1, 2, 3, 1)], [("Fab", "HC:1", "LC1:1", None), ("Fab", "HC:2", "LC2:1", None)]),
        (DVD_IG, [(1, 1, 2, 1), (1, 2, 2, 2)], [("Fv", "HC:1", "LC:1", None), ("Fab", "HC:2", "LC:2", None)]),
        # pairs are not changed by the format, e.g. crossed V domains of CODV-Ig
        (DVD_IG, [(1, 1, 2, 2), (1, 2, 2, 1)], [("Fv", "HC:1", "LC:2", None), ("Fv", "HC:2", "LC:1", None)]),
        (
            {"LC-HC": VK + CK + VH2 + IGHG1, "Fd": VH + CH1, "LC": VK2 + CK},
            [(2, 1, 1, 1), (1, 2, 3, 1)],
            [("Fab", "Fd:1", "LC-HC:1", None), ("Fab", "LC-HC:2", "LC:1", None)],
        ),
        (
            {"HC": VH + IGHG1 + LINKER + VH2 + LINKER + VK2, "LC": VK + CK},
            [(1, 1, 2, 1), (1, 2, 1, 3)],
            [("Fab", "HC:1", "LC:1", None), ("scFv", "HC:2", "HC:3", 15)],
        ),
        (
            {"BiTE": VK + LINKER + VH + "GGGGS" + VH2 + LINKER + VK2},
            [(1, 2, 1, 1), (1, 3, 1, 4)],
            [("scFv", "BiTE:2", "BiTE:1", 15), ("scFv", "BiTE:3", "BiTE:4", 15)],
        ),
        ({"chain": VH + "GGGGS" + VK}, [(1, 1, 1, 2)], [("diabody", "chain:1", "chain:2", 5)]),
        (
            {"A": VH + "GGGGS" + VK2, "B": VH2 + "GGGGS" + VK},
            [(1, 1, 2, 2), (2, 1, 1, 2)],
            [("diabody", "A:1", "B:2", None), ("diabody", "B:1", "A:2", None)],
        ),
        (
            {"tandem": VH + "GGGGS" + VK2 + "GGGGS" + VH2 + "GGGGS" + VK},
            [(1, 1, 1, 4), (1, 3, 1, 2)],
            [("diabody", "tandem:1", "tandem:4", None), ("diabody", "tandem:3", "tandem:2", 5)],
        ),
        ({"scFab": VK + CK + 6 * "GGGGS" + VH + CH1}, [(1, 2, 1, 1)], [("Fab", "scFab:2", "scFab:1", None)]),
        ({"H": VH, "L": VK}, [(1, 1, 2, 1)], [("Fv", "H:1", "L:1", None)]),
        ({"HC": VH + IGHG1, "L": VK}, [(1, 1, 2, 1)], [("Fv", "HC:1", "L:1", None)]),
        # V domains of a scFv paired with other V domains
        (
            {"A": VH + LINKER + VK2, "B": VK},
            [(1, 1, 2, 1)],
            [("Fv", "A:1", "B:1", None), ("unpaired", None, "A:2", None)],
        ),
        # V domains not given in pairs are single domains or unpaired
        (
            {"HC1": VH + IGHG1, "HC2": VH2 + IGHG1, "LC": VK + CK, "VHH-Fc": VHH + FC},
            [(1, 1, 3, 1)],
            [
                ("Fab", "HC1:1", "LC:1", None),
                ("unpaired", "HC2:1", None, None),
                ("single domain", "VHH-Fc:1", None, None),
            ],
        ),
        ({"HC": VH + IGHG1, "LC": VK + CK}, [], [("unpaired", "HC:1", None, None), ("unpaired", None, "LC:1", None)]),
    ],
    ids=[
        "IgG",
        "bispecific",
        "chain-IDs",
        "common-light-chain",
        "common-heavy-chain",
        "CrossMab",
        "tandem-Fab",
        "DVD-Ig",
        "DVD-Ig-crossed",
        "FIT-Ig",
        "IgG-scFv",
        "BiTE",
        "diabody",
        "heterodimeric-diabody",
        "tandem-diabody",
        "scFab",
        "V-domains",
        "HC-with-VL",
        "scFv-V-domains",
        "not-given",
        "none",
    ],
)
def test_pairs(riot_aa: RiotNumberingAA, chains: dict[str, str], pairs: list, expected: list):
    molecule = riot_aa.annotate_molecule("molecule", chains, pairs=pairs)
    assert numbered_binding_sites(molecule) == expected
    assert not any(binding_site.ambiguous or binding_site.candidates for binding_site in molecule.binding_sites)


@pytest.mark.parametrize(
    "chains, pairs, error",
    [
        (BISPECIFIC, [(5, 1, 3, 1)], r"No chain 5 in V domain pair \(5, 1, 3, 1\), the molecule has 4 chains"),
        (BISPECIFIC, [(0, 1, 3, 1)], "No chain 0"),
        (BISPECIFIC, [("HC3", 1, "LC1", 1)], "Unknown chain HC3"),
        (
            BISPECIFIC,
            [(1, 2, 3, 1)],
            r"No V domain 2 in chain 1 \(HC1\) in V domain pair \(1, 2, 3, 1\), "
            r"V domains of the chain: 1 VH \(layout VH-CH1-H-CH2-CH3\)",
        ),
        (
            {"HC": VH + IGHG1, "Fc": FC, "LC": VK + CK},
            [(1, 1, 2, 1)],
            r"No V domain 1 in chain 2 \(Fc\) .* V domains of the chain: none \(layout H-CH2-CH3\)",
        ),
        (
            BISPECIFIC,
            [(1, 1, 2, 1)],
            r"pairs two heavy V domains: VH 1 of chain 1 \(HC1\) and VH 1 of chain 2 \(HC2\)",
        ),
        (BISPECIFIC, [(3, 1, 4, 1)], "pairs two light V domains"),
        (BISPECIFIC, [(1, 1, 1, 1)], r"pairs VH 1 of chain 1 \(HC1\) with itself"),
        (BISPECIFIC, [(1, 1, 3, 1), (3, 1, 1, 1)], r"Duplicate V domain pair \(3, 1, 1, 1\)"),
        (BISPECIFIC, [(1, 1, 3)], r"is not \(chain, V domain number, chain, V domain number\)"),
        (BISPECIFIC, [(1, "1", 3, 1)], r"is not \(chain, V domain number, chain, V domain number\)"),
    ],
    ids=[
        "chain-number",
        "chain-number-0",
        "chain-ID",
        "V-domain-number",
        "chain-without-V-domains",
        "two-heavy",
        "two-light",
        "same-V-domain",
        "duplicate",
        "short",
        "not-number",
    ],
)
def test_invalid_pairs(riot_aa: RiotNumberingAA, chains: dict[str, str], pairs: list, error: str):
    with pytest.raises(ValueError, match=error):
        riot_aa.annotate_molecule("molecule", chains, pairs=pairs)


@pytest.mark.parametrize(
    "chains, pairs, expected",
    [
        # common light chain is in binding sites of both heavy chains
        (
            {"HC1": VH + IGHG1, "HC2": VH2 + IGHG1, "LC": VK + CK},
            None,
            {
                "HC1:1": (1, "Fab", "LC", 1, False, None),
                "HC2:1": (2, "Fab", "LC", 1, False, None),
                "LC:1": ("1;2", "Fab;Fab", "HC1;HC2", "1;1", False, None),
            },
        ),
        # light chains are candidates of ambiguous binding sites of heavy chains
        (
            BISPECIFIC,
            None,
            {
                "HC1:1": (1, "Fab", None, None, True, "LC1:1;LC2:1"),
                "HC2:1": (2, "Fab", None, None, True, "LC1:1;LC2:1"),
                "LC1:1": ("1;2", "Fab;Fab", None, None, True, "HC1:1;HC2:1"),
                "LC2:1": ("1;2", "Fab;Fab", None, None, True, "HC1:1;HC2:1"),
            },
        ),
        # common heavy chain is in binding sites of both light chains
        (
            {"HC": VH + IGHG1, "LC1": VK + CK, "LC2": VL + C_GERMLINES_AA["IGLC2"]},
            [(1, 1, 2, 1), (1, 1, 3, 1)],
            {
                "HC:1": ("1;2", "Fab;Fab", "LC1;LC2", "1;1", False, None),
                "LC1:1": (1, "Fab", "HC", 1, False, None),
                "LC2:1": (2, "Fab", "HC", 1, False, None),
            },
        ),
    ],
    ids=["common-light-chain", "bispecific", "common-heavy-chain"],
)
def test_layout_binding_sites(riot_aa: RiotNumberingAA, chains: dict[str, str], pairs: Optional[list], expected: dict):
    molecule = riot_aa.annotate_molecule("molecule", chains, pairs=pairs)
    v_domain_rows = {
        f"{row['chain_id']}:{row['domain_index']}": (
            row["binding_site_id"],
            row["binding_site_type"],
            row["partner_chain_id"],
            row["partner_domain_index"],
            row["ambiguous"],
            row["candidates"],
        )
        for row in molecule_layout_rows(molecule)
        if row["rearrangement_index"] is not None and row["domain"] in ("VH", "VK", "VL")
    }
    assert v_domain_rows == expected


def test_chain_list(riot_aa: RiotNumberingAA):
    # chains are numbered, chain numbers are chain IDs
    molecule = riot_aa.annotate_molecule("molecule", [VH + IGHG1, VK + CK], pairs=[(1, 1, 2, 1)])
    assert [chain.sequence_header for chain in molecule.chains.values()] == ["molecule|1", "molecule|2"]
    assert numbered_binding_sites(molecule) == [("Fab", "1:1", "2:1", None)]

    with pytest.raises(TypeError):
        riot_aa.annotate_molecule("molecule", VH + IGHG1)


def test_v_domain_numbers(riot_aa: RiotNumberingAA):
    molecule = riot_aa.annotate_molecule("molecule", {"HC": VH + IGHG1 + LINKER + VH2 + LINKER + VK2, "LC": VK + CK})

    assert [
        (row["chain_id"], row["domain"], row["domain_index"], row["v_domain_number"], row["partner_v_domain_number"])
        for row in molecule_layout_rows(molecule)
        if row["v_domain_number"]
    ] == [("HC", "VH", 1, 1, 1), ("HC", "VH", 6, 2, 3), ("HC", "VK", 7, 3, 2), ("LC", "VK", 1, 1, 1)]
    numbers = [(row["heavy_v_domain_number"], row["light_v_domain_number"]) for row in binding_site_rows(molecule)]
    assert numbers == [(1, 1), (2, 3)]


@pytest.mark.parametrize("reverse", [False, True])
def test_nucleotide_pairs(reverse: bool):
    riot_nt = get_or_create_riot_nt(return_all_domains=True)
    vh = nt_rearrangement("IGHV3-23*01", "GGCTAC", "IGHJ4*02", 2)
    vk = nt_rearrangement("IGKV1-27*01", "CAGCAGAGTTACAGTACCCCT", "IGKJ1*01", 1)
    scfv = vh + 3 * "GGTGGCGGTGGCTCG" + vk
    molecule = riot_nt.annotate_molecule(
        "molecule",
        [
            nt_rearrangement("IGHV3-23*01", "GGCTAC", "IGHJ4*02", 2, "IGHG1"),
            nt_rearrangement("IGKV1-27*01", "CAGCAGAGTTACAGTACCCCT", "IGKJ1*01", 1, "IGKC"),
            reverse_complement(scfv) if reverse else scfv,
        ],
        pairs=[(1, 1, 2, 1), (3, 1, 3, 2)],
    )

    # V domains are numbered in the reverse complemented sequence of chains on the reverse strand
    assert molecule.chains["3"].rev_comp == reverse
    assert numbered_binding_sites(molecule) == [("Fab", "1:1", "2:1", None), ("scFv", "3:1", "3:2", 15)]


def test_read_v_domain_pairs(tmp_path: Path):
    pairs_path = tmp_path / "pairs.tsv"
    pairs_path.write_text(
        "molecule_id\tchain_1\tv_domain_1\tchain_2\tv_domain_2\n"
        "bispecific\t1\t1\t2\t1\n"
        "IgG\t1\t1\t2\t1\n"
        "bispecific\t3\t1\t4\t1\n"
    )
    assert read_v_domain_pairs(pairs_path) == {"bispecific": [(1, 1, 2, 1), (3, 1, 4, 1)], "IgG": [(1, 1, 2, 1)]}

    pairs_path = tmp_path / "pairs.csv"
    pairs_path.write_text("molecule_id,chain_1,v_domain_1,chain_2\nIgG,1,1,2\n")
    with pytest.raises(ValueError, match="Missing columns .*: v_domain_2"):
        read_v_domain_pairs(pairs_path)
    pairs_path.write_text("molecule_id,chain_1,v_domain_1,chain_2,v_domain_2\nIgG,1,1,2,1\nIgG2,HC,1,LC,1\n")
    with pytest.raises(ValueError, match="Invalid V domain pair in line 3"):
        read_v_domain_pairs(pairs_path)


def test_nucleotide_molecule():
    molecule = get_or_create_riot_nt(return_all_domains=True).annotate_molecule(
        "molecule",
        {
            "HC": nt_rearrangement("IGHV3-23*01", "GGCTAC", "IGHJ4*02", 2, "IGHG1"),
            "LC": nt_rearrangement("IGKV1-27*01", "CAGCAGAGTTACAGTACCCCT", "IGKJ1*01", 1, "IGKC"),
        },
    )

    assert {chain_id: chain.layout for chain_id, chain in molecule.chains.items()} == {
        "HC": "VH-CH1-H-CH2-CH3",
        "LC": "VK-CL",
    }
    assert binding_sites(molecule) == [("Fab", "HC:1", "LC:1", None)]


# J region of VH mutated (YFDYWGQGTLVTVSS -> YFAAWGAGTAVTASS), J is not found, the J motif (W-G-x-G) is kept
VH_WITHOUT_J = VH.replace("YFDYWGQGTLVTVSS", "YFAAWGAGTAVTASS")


@pytest.mark.parametrize(
    "chains, expected, linker_length",
    [
        ({"chain": VH_WITHOUT_J + "GGGGS" + VK}, ("diabody", "chain:1", "chain:2", None), 5),
        ({"chain": VH_WITHOUT_J + LINKER + VK}, ("scFv", "chain:1", "chain:2", None), len(LINKER)),
    ],
    ids=["diabody", "scFv"],
)
def test_v_domains_without_j(riot_aa: RiotNumberingAA, chains: dict[str, str], expected: tuple, linker_length: int):
    molecule = riot_aa.annotate_molecule("molecule", chains)

    assert molecule.chains["chain"].rearrangements[0].j_call is None
    assert binding_sites(molecule) == [expected]
    assert molecule.binding_sites[0].linker_length == linker_length


@pytest.mark.parametrize(
    "linker, expected_type", [("GGTGGCGGTGGCTCG", "diabody"), (2 * "GGTGGCGGTGGCTCG" + "GGTGGC", "scFv")]
)
def test_nucleotide_linker_length(linker: str, expected_type: str):
    # J germlines end with the nucleotide of the first C codon, which matches the first nucleotide of the linker
    vh = nt_rearrangement("IGHV3-23*01", "GGCTAC", "IGHJ4*02", 2)
    vk = nt_rearrangement("IGKV1-27*01", "CAGCAGAGTTACAGTACCCCT", "IGKJ1*01", 1)
    molecule = get_or_create_riot_nt(return_all_domains=True).annotate_molecule("molecule", {"chain": vh + linker + vk})

    assert [
        (binding_site.binding_site_type, binding_site.linker_length) for binding_site in molecule.binding_sites
    ] == [(expected_type, len(linker) // 3)]


def test_molecule_output(riot_aa: RiotNumberingAA, tmp_path: Path):
    molecules = [
        riot_aa.annotate_molecule("IgG", {"HC": VH + IGHG1, "LC": VK + CK}, c_scheme=CScheme.EU),
        riot_aa.annotate_molecule("VHH-Fc", {"chain": VHH + FC}),
    ]
    write_molecule_layout_csv(tmp_path / "layout.csv", molecules)
    write_binding_sites_csv(tmp_path / "binding_sites.csv", molecules)
    write_molecules_json(tmp_path / "molecules.json", molecules)

    with (tmp_path / "layout.csv").open() as layout_file:
        layout = list(csv.DictReader(layout_file))
    assert [(row["molecule_id"], row["chain_id"], row["domain"]) for row in layout] == [
        ("IgG", "HC", "VH"),
        ("IgG", "HC", "CH1"),
        ("IgG", "HC", "H"),
        ("IgG", "HC", "CH2"),
        ("IgG", "HC", "CH3"),
        ("IgG", "LC", "VK"),
        ("IgG", "LC", "CL"),
        ("VHH-Fc", "chain", "VH"),
        ("VHH-Fc", "chain", "H"),
        ("VHH-Fc", "chain", "CH2"),
        ("VHH-Fc", "chain", "CH3"),
    ]
    vh, ch1 = layout[0], layout[1]
    assert (vh["locus"], vh["binding_site_type"], vh["partner_chain_id"], vh["partner_domain_index"]) == (
        "igh",
        "Fab",
        "LC",
        "1",
    )
    assert (ch1["c_call"], ch1["rearrangement_index"], ch1["binding_site_id"]) == ("IGHG1", "0", "")

    with (tmp_path / "binding_sites.csv").open() as binding_sites_file:
        rows = list(csv.DictReader(binding_sites_file))
    assert [(row["molecule_id"], row["binding_site_type"], row["light_locus"]) for row in rows] == [
        ("IgG", "Fab", "igk"),
        ("VHH-Fc", "single domain", ""),
    ]

    with (tmp_path / "molecules.json").open() as json_file:
        molecules_json = json.load(json_file)
    assert [molecule["molecule_id"] for molecule in molecules_json] == ["IgG", "VHH-Fc"]
    assert molecules_json[0]["chains"]["HC"]["constant_domains"][2]["c_scheme_residue_mapping"]["297"] == "N"


def fab_arm(chain_id: str, heavy: bool) -> _FabArm:
    v_domain = _VDomain(0, chain_id, 0, ChainDomain("VH" if heavy else "VK", 1, 120, rearrangement_index=0))
    return _FabArm((v_domain,), ChainDomain("CH1" if heavy else "CL", 121, 220), "CH1" if heavy else "CL")


def test_one_to_one_partners():
    ch1_arms = [fab_arm(f"HC{i}", True) for i in range(10)]
    cl_arms = [fab_arm(f"LC{i}", False) for i in range(10)]

    # every Fab arm can pair with every other, so all partners are ambiguous
    all_partners = {ch1_arm: list(cl_arms) for ch1_arm in ch1_arms}
    assert _one_to_one_partners(ch1_arms, cl_arms, all_partners) == all_partners

    # FIT-Ig like: HC1 pairs with LC1 only, so HC0 pairs with LC0
    partners = {ch1_arms[0]: cl_arms[:2], ch1_arms[1]: cl_arms[1:2]}
    assert _one_to_one_partners(ch1_arms[:2], cl_arms[:2], partners) == {
        ch1_arms[0]: cl_arms[:1],
        ch1_arms[1]: cl_arms[1:2],
    }

    # no one-to-one pairing, the last CL arm has no partner
    no_pairing = {ch1_arm: cl_arms[:-1] for ch1_arm in ch1_arms}
    assert _one_to_one_partners(ch1_arms, cl_arms, no_pairing) is None
    assert _one_to_one_partners(ch1_arms, cl_arms[:-1], no_pairing) is None
