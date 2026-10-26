import csv
import json
from pathlib import Path
from typing import Optional

import pytest

from riot_na import CScheme, MoleculeAnnotation
from riot_na.api.molecule_annotation import (
    _FabArm,
    _one_to_one_partners,
    _VDomain,
    binding_site_rows,
)
from riot_na.api.riot_numbering import (
    RiotNumberingAA,
    get_or_create_riot_aa,
    get_or_create_riot_nt,
)
from riot_na.common.io import (
    write_binding_sites_csv,
    write_molecule_layout_csv,
    write_molecules_json,
)
from riot_na.data.model import ChainDomain
from tests.test_c_numbering import (
    C_GERMLINES_AA,
    VH,
    VK,
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


@pytest.mark.parametrize(
    "chains, pairs, expected",
    [
        (
            BISPECIFIC,
            [("HC1", "LC2"), ("HC2", "LC1")],
            [("Fab", "HC1:1", "LC2:1", None), ("Fab", "HC2:1", "LC1:1", None)],
        ),
        # HC2 pairs with the light chain left
        (BISPECIFIC, [("HC1", "LC2")], [("Fab", "HC1:1", "LC2:1", None), ("Fab", "HC2:1", "LC1:1", None)]),
        # light chains given in pairs are not reserved, they can be common light chains
        (
            {"HC1": VH + IGHG1, "HC2": VH2 + IGHG1, "LC": VK + CK},
            [("HC1", "LC")],
            [("Fab", "HC1:1", "LC:1", None), ("Fab", "HC2:1", "LC:1", None)],
        ),
        (
            TANDEM_FAB,
            None,
            [("Fab", "HC:1", None, "LC1:1;LC2:1"), ("Fab", "HC:3", None, "LC1:1;LC2:1")],
        ),
        # Fab arms of one chain given as V domains (chain ID and domain index)
        (
            TANDEM_FAB,
            [(("HC", 1), "LC1"), (("HC", 3), "LC2")],
            [("Fab", "HC:1", "LC1:1", None), ("Fab", "HC:3", "LC2:1", None)],
        ),
        # V domain only chains
        ({"H": VH, "L": VK}, None, [("single domain", "H:1", None, None), ("unpaired", None, "L:1", None)]),
        ({"H": VH, "L": VK}, [("H", "L")], [("Fv", "H:1", "L:1", None)]),
        # heavy chain without its light chain constant domain
        ({"HC": VH + IGHG1, "L": VK}, [("HC", "L")], [("Fv", "HC:1", "L:1", None)]),
        (
            DIABODIES,
            None,
            [
                ("diabody", "A:1", "A:2", "B:1;B:2;C:1;C:2"),
                ("diabody", "B:1", "B:2", "A:1;A:2;C:1;C:2"),
                ("diabody", "C:1", "C:2", "A:1;A:2;B:1;B:2"),
            ],
        ),
        (
            DIABODIES,
            [("A", "B")],
            [("diabody", "A:1", "B:2", None), ("diabody", "B:1", "A:2", None), ("diabody", "C:1", "C:2", None)],
        ),
    ],
    ids=[
        "bispecific",
        "one-pair",
        "common-light-chain",
        "tandem-Fab",
        "tandem-Fab-domains",
        "V-domains",
        "Fv",
        "HC-with-VL",
        "diabodies",
        "heterodimeric-diabody",
    ],
)
def test_pairs(riot_aa: RiotNumberingAA, chains: dict[str, str], pairs: Optional[list], expected: list):
    molecule = riot_aa.annotate_molecule("molecule", chains, pairs=pairs)
    assert binding_sites(molecule) == expected


def test_invalid_pairs(riot_aa: RiotNumberingAA):
    with pytest.raises(ValueError, match="Unknown chain"):
        riot_aa.annotate_molecule("molecule", BISPECIFIC, pairs=[("HC1", "LC3")])
    # domain 2 of HC is CH1
    with pytest.raises(ValueError, match="No V domain with domain index 2 in chain HC"):
        riot_aa.annotate_molecule("molecule", TANDEM_FAB, pairs=[(("HC", 2), "LC1")])


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
