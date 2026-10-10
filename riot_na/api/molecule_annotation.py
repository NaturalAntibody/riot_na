"""
Molecule annotation: binding sites formed by V domains of the chains of a molecule.

V domains are paired by the format of the chains, in this order:
- scFv: adjacent heavy and light V domains of a chain (no constant domain in between) joined by a linker of at least
  MIN_SCFV_LINKER_LENGTH residues,
- diabody: the same with a shorter linker, too short to pair V domains of one chain. V domains pair with V domains of
  another copy of the chain (homodimer) or, if two chains have such V domains, of the other chain (heterodimer),
- Fab: V domain followed by CH1 pairs with V domain followed by CL, if one of them is heavy and the other light
  (which pairs domain crossovers of CrossMabs too). V domains stacked in front of them (e.g. outer V domains
  of DVD-Ig) pair by their position. Fab arms of one chain pair only if joined by a linker of at least
  MIN_SCFAB_LINKER_LENGTH residues (single chain Fab). Fab arms with more than one possible partner are ambiguous,
  unless the chains to pair are given,
- single domain: other heavy V domains (e.g. VHH), other light V domains are unpaired.
"""

from dataclasses import dataclass, field
from typing import Any, Mapping, Optional, Sequence

from riot_na.data.model import (
    BindingSite,
    ChainAnnotation,
    ChainDomain,
    MoleculeAnnotation,
    VDomainRef,
)

MIN_SCFV_LINKER_LENGTH = 12
MIN_SCFAB_LINKER_LENGTH = 20
# shorter CH1 or CL fragments (e.g. CH1 derived linkers of DVD-Ig) don't form Fab arms, unless at chain end
MIN_FAB_CONSTANT_DOMAIN_LENGTH = 30
HEAVY_V_DOMAIN = "VH"
V_DOMAINS = (HEAVY_V_DOMAIN, "VK", "VL")
CH1 = "CH1"
CL = "CL"
FAB_CONSTANT_DOMAINS = (CH1, CL)
# limit of one-to-one pairings of Fab arms enumerated to find unambiguous partners
MAX_ONE_TO_ONE_PAIRINGS = 1000

FAB = "Fab"
STACKED_FV = "Fv"
SCFV = "scFv"
DIABODY = "diabody"
SINGLE_DOMAIN = "single domain"
UNPAIRED = "unpaired"

LAYOUT_FIELDS = (
    "molecule_id",
    "chain_id",
    "chain_layout",
    "domain_index",
    "domain",
    "sequence_start",
    "sequence_end",
    "locus",
    "species",
    "v_call",
    "j_call",
    "c_call",
    "c_identity",
    "rearrangement_index",
    "binding_site_id",
    "binding_site_type",
    "partner_chain_id",
    "partner_domain_index",
    "ambiguous",
    "candidates",
)
BINDING_SITE_FIELDS = (
    "molecule_id",
    "binding_site_id",
    "binding_site_type",
    "heavy_chain_id",
    "heavy_domain_index",
    "heavy_v_call",
    "light_chain_id",
    "light_domain_index",
    "light_locus",
    "light_v_call",
    "linker_length",
    "ambiguous",
    "candidates",
)


@dataclass(frozen=True)
class _VDomain:
    chain_order: int
    chain_id: str
    # position in chain domains
    position: int
    domain: ChainDomain = field(compare=False)

    @property
    def heavy(self) -> bool:
        return self.domain.domain == HEAVY_V_DOMAIN

    @property
    def ref(self) -> VDomainRef:
        assert self.domain.rearrangement_index is not None
        return VDomainRef(self.chain_id, self.domain.rearrangement_index)


@dataclass(frozen=True)
class _FabArm:
    # V domains stacked in front of the constant domain, the last one is followed by it
    v_domains: tuple[_VDomain, ...]
    constant_domain: ChainDomain = field(compare=False)

    @property
    def chain_id(self) -> str:
        return self.v_domains[-1].chain_id

    @property
    def heavy(self) -> bool:
        return self.v_domains[-1].heavy


class _BindingSitePairing:  # pylint: disable=too-few-public-methods
    def __init__(
        self, chains: Mapping[str, ChainAnnotation], residue_length: int, pairs: Optional[Sequence[tuple[str, str]]]
    ) -> None:
        self.chains = chains
        self.residue_length = residue_length
        self.partner_chains: dict[str, set[str]] = {}
        for chain_id, partner_chain_id in pairs or ():
            for chain in (chain_id, partner_chain_id):
                if chain not in chains:
                    raise ValueError(f"Unknown chain in pairs: {chain}")
            self.partner_chains.setdefault(chain_id, set()).add(partner_chain_id)
            self.partner_chains.setdefault(partner_chain_id, set()).add(chain_id)

        self.v_domains = {
            chain_id: [
                _VDomain(chain_order, chain_id, position, domain)
                for position, domain in enumerate(chain.domains)
                if domain.domain in V_DOMAINS
            ]
            for chain_order, (chain_id, chain) in enumerate(chains.items())
        }
        self.binding_sites: list[BindingSite] = []
        self.paired: set[_VDomain] = set()

    def pair(self) -> list[BindingSite]:
        self._pair_single_chain_fvs()
        self._pair_fab_arms()

        for chain_v_domains in self.v_domains.values():
            for v_domain in chain_v_domains:
                if v_domain not in self.paired:
                    self._add_binding_site(SINGLE_DOMAIN if v_domain.heavy else UNPAIRED, [v_domain])

        domain_order = {
            v_domain.ref: (v_domain.chain_order, v_domain.position)
            for chain_v_domains in self.v_domains.values()
            for v_domain in chain_v_domains
        }
        binding_sites = sorted(
            self.binding_sites,
            key=lambda binding_site: min(
                domain_order[domain] for domain in (binding_site.heavy_domain, binding_site.light_domain) if domain
            ),
        )
        for binding_site_id, binding_site in enumerate(binding_sites, start=1):
            binding_site.binding_site_id = binding_site_id
        return binding_sites

    def _linker_length(self, first: ChainDomain, second: ChainDomain) -> int:
        return (second.sequence_start - first.sequence_end - 1) // self.residue_length

    def _domain_length(self, domain: ChainDomain) -> int:
        return (domain.sequence_end - domain.sequence_start + 1) // self.residue_length

    def _add_binding_site(  # pylint: disable=too-many-arguments
        self,
        binding_site_type: str,
        v_domains: Sequence[_VDomain],
        ambiguous: bool = False,
        candidates: Sequence[_VDomain] = (),
        linker_length: Optional[int] = None,
    ) -> None:
        heavy_domains = [v_domain.ref for v_domain in v_domains if v_domain.heavy]
        light_domains = [v_domain.ref for v_domain in v_domains if not v_domain.heavy]
        self.binding_sites.append(
            BindingSite(
                binding_site_id=0,
                binding_site_type=binding_site_type,
                heavy_domain=heavy_domains[0] if heavy_domains else None,
                light_domain=light_domains[0] if light_domains else None,
                ambiguous=ambiguous,
                candidates=[candidate.ref for candidate in candidates],
                linker_length=linker_length,
            )
        )
        self.paired.update(v_domains)

    def _pair_single_chain_fvs(self) -> None:
        """Pair adjacent heavy and light V domains of chains: scFv or, joined by short linkers, diabodies."""
        diabody_v_domains: dict[str, list[tuple[_VDomain, _VDomain, int]]] = {}

        for chain_id, chain_v_domains in self.v_domains.items():
            i = 0
            while i + 1 < len(chain_v_domains):
                first, second = chain_v_domains[i], chain_v_domains[i + 1]
                if second.position != first.position + 1 or first.heavy == second.heavy:
                    i += 1
                    continue

                linker_length = self._linker_length(first.domain, second.domain)
                if linker_length >= MIN_SCFV_LINKER_LENGTH:
                    self._add_binding_site(SCFV, [first, second], linker_length=linker_length)
                else:
                    diabody_v_domains.setdefault(chain_id, []).append((first, second, linker_length))
                    self.paired.update((first, second))
                i += 2

        heterodimer = len(diabody_v_domains) == 2 and all(len(pairs) == 1 for pairs in diabody_v_domains.values())
        if heterodimer:
            (first, second, _), (other_first, other_second, _) = [pairs[0] for pairs in diabody_v_domains.values()]
            # every V domain pairs with the complementary V domain of the other chain
            self._add_binding_site(DIABODY, [first, other_second if first.heavy != other_second.heavy else other_first])
            self._add_binding_site(
                DIABODY, [second, other_first if second.heavy != other_first.heavy else other_second]
            )
            return

        for chain_id, pairs in diabody_v_domains.items():
            for first, second, linker_length in pairs:
                # V domains of other chains could pair with them too
                candidates = [
                    v_domain
                    for other_chain_id, other_pairs in diabody_v_domains.items()
                    if other_chain_id != chain_id
                    for other_pair in other_pairs
                    for v_domain in other_pair[:2]
                ]
                self._add_binding_site(
                    DIABODY,
                    [first, second],
                    ambiguous=bool(candidates),
                    candidates=candidates,
                    linker_length=linker_length,
                )

    def _fab_arms(self) -> list[_FabArm]:
        """V domains followed by CH1 or CL, with V domains stacked in front of them."""
        fab_arms = []
        for chain_id, chain in self.chains.items():
            domains = chain.domains
            v_domains = {v_domain.position: v_domain for v_domain in self.v_domains[chain_id]}

            for position, domain in enumerate(domains):
                v_domain = v_domains.get(position - 1)
                full_constant_domain = (
                    self._domain_length(domain) >= MIN_FAB_CONSTANT_DOMAIN_LENGTH or position == len(domains) - 1
                )
                if domain.domain not in FAB_CONSTANT_DOMAINS or not full_constant_domain:
                    continue
                if v_domain is None or v_domain in self.paired:
                    continue

                stack = [v_domain]
                stack_position = position - 2
                while stack_position >= 0:
                    stacked_domain = domains[stack_position]
                    stacked_v_domain = v_domains.get(stack_position)
                    if stacked_v_domain is None:
                        # short constant domain fragments can join stacked V domains
                        if self._domain_length(stacked_domain) >= MIN_FAB_CONSTANT_DOMAIN_LENGTH:
                            break
                    elif stacked_v_domain.heavy != v_domain.heavy or stacked_v_domain in self.paired:
                        break
                    else:
                        stack.insert(0, stacked_v_domain)
                    stack_position -= 1

                fab_arms.append(_FabArm(tuple(stack), domain))
        return fab_arms

    def _fab_partners(self, fab_arm: _FabArm, fab_arms: Sequence[_FabArm]) -> list[_FabArm]:
        partners = []
        for other_arm in fab_arms:
            if (
                other_arm.constant_domain.domain == fab_arm.constant_domain.domain
                or other_arm.heavy == fab_arm.heavy
                or len(other_arm.v_domains) != len(fab_arm.v_domains)
            ):
                continue
            if (
                fab_arm.chain_id in self.partner_chains
                and other_arm.chain_id not in self.partner_chains[fab_arm.chain_id]
            ):
                continue
            if other_arm.chain_id == fab_arm.chain_id:
                first, second = sorted((fab_arm, other_arm), key=lambda arm: arm.constant_domain.sequence_start)
                # single chain Fab
                if self._linker_length(first.constant_domain, second.v_domains[0].domain) < MIN_SCFAB_LINKER_LENGTH:
                    continue
            partners.append(other_arm)
        return partners

    def _pair_fab_arms(self) -> None:
        fab_arms = self._fab_arms()
        ch1_arms = [fab_arm for fab_arm in fab_arms if fab_arm.constant_domain.domain == CH1]
        cl_arms = [fab_arm for fab_arm in fab_arms if fab_arm.constant_domain.domain == CL]
        fab_arm_partners = {fab_arm: self._fab_partners(fab_arm, cl_arms) for fab_arm in ch1_arms}
        fab_arm_partners = _one_to_one_partners(ch1_arms, cl_arms, fab_arm_partners) or fab_arm_partners
        cl_arms_with_partners: set[_FabArm] = set()

        for fab_arm in ch1_arms:
            partners = fab_arm_partners[fab_arm]
            cl_arms_with_partners.update(partners)
            stack_size = len(fab_arm.v_domains)

            for position, v_domain in enumerate(fab_arm.v_domains):
                binding_site_type = FAB if position == stack_size - 1 else STACKED_FV
                if len(partners) == 1:
                    self._add_binding_site(binding_site_type, [v_domain, partners[0].v_domains[position]])
                elif partners:
                    self._add_binding_site(
                        binding_site_type,
                        [v_domain],
                        ambiguous=True,
                        candidates=[partner.v_domains[position] for partner in partners],
                    )
                else:
                    self._add_binding_site(UNPAIRED, [v_domain])

        for fab_arm in cl_arms:
            if fab_arm not in cl_arms_with_partners:
                for v_domain in fab_arm.v_domains:
                    if v_domain not in self.paired:
                        self._add_binding_site(UNPAIRED, [v_domain])
        # V domains of CL arms with ambiguous partners are listed as candidates of the CH1 arms
        for fab_arm in cl_arms_with_partners:
            self.paired.update(fab_arm.v_domains)


def _one_to_one_partners(
    ch1_arms: Sequence[_FabArm], cl_arms: Sequence[_FabArm], fab_arm_partners: Mapping[_FabArm, list[_FabArm]]
) -> Optional[dict[_FabArm, list[_FabArm]]]:
    """
    If there are as many CH1 as CL Fab arms, every CL arm likely pairs with one CH1 arm (e.g. FIT-Ig). Partners
    of CH1 arms are then the CL arms they are paired with in any one-to-one pairing. Returns None, if the numbers
    of arms differ (e.g. common light chain pairs with many heavy chains) or there is no one-to-one pairing.
    """
    if len(ch1_arms) != len(cl_arms):
        return None

    pairings: list[dict[_FabArm, _FabArm]] = []

    def pair_from(arm_index: int, pairing: dict[_FabArm, _FabArm]):
        if len(pairings) >= MAX_ONE_TO_ONE_PAIRINGS:
            return
        if arm_index == len(ch1_arms):
            pairings.append(dict(pairing))
            return
        fab_arm = ch1_arms[arm_index]
        for partner in fab_arm_partners[fab_arm]:
            if partner not in pairing.values():
                pairing[fab_arm] = partner
                pair_from(arm_index + 1, pairing)
                del pairing[fab_arm]

    pair_from(0, {})
    if not pairings:
        return None
    return {
        fab_arm: [
            partner for partner in fab_arm_partners[fab_arm] if any(pairing[fab_arm] == partner for pairing in pairings)
        ]
        for fab_arm in ch1_arms
    }


def pair_binding_sites(
    chains: Mapping[str, ChainAnnotation],
    residue_length: int = 1,
    pairs: Optional[Sequence[tuple[str, str]]] = None,
) -> list[BindingSite]:
    """
    Binding sites of the V domains of a molecule chains.

    - residue_length: length of residues in chain sequence positions, 3 for nucleotide chains.
    - pairs: chains which Fab arms pair, other chains are not considered as their partners. They resolve ambiguous
    Fab arms, e.g. of bispecific antibodies with two different heavy and light chains.
    """
    return _BindingSitePairing(chains, residue_length, pairs).pair()


def molecule_layout_rows(molecule: MoleculeAnnotation) -> list[dict[str, Any]]:
    """Domains of all chains of a molecule, one row per domain, with binding sites of V domains."""
    domain_indices = _domain_indices(molecule)
    binding_sites: dict[VDomainRef, BindingSite] = {}
    for binding_site in molecule.binding_sites:
        for v_domain in (binding_site.heavy_domain, binding_site.light_domain):
            if v_domain:
                binding_sites[v_domain] = binding_site

    rows = []
    for chain_id, chain in molecule.chains.items():
        for domain_index, domain in enumerate(chain.domains, start=1):
            row: dict[str, Any] = {
                "molecule_id": molecule.molecule_id,
                "chain_id": chain_id,
                "chain_layout": chain.layout,
                "domain_index": domain_index,
                "domain": domain.domain,
                "sequence_start": domain.sequence_start,
                "sequence_end": domain.sequence_end,
                "rearrangement_index": domain.rearrangement_index,
            }

            if domain.constant_domain_index is not None:
                constant_domain = chain.constant_domains[domain.constant_domain_index]
                row |= {
                    "locus": constant_domain.locus,
                    "species": constant_domain.species,
                    "c_call": constant_domain.c_call,
                    "c_identity": constant_domain.identity_aa,
                }
            else:
                assert domain.rearrangement_index is not None
                rearrangement = chain.rearrangements[domain.rearrangement_index]
                row |= {
                    "locus": rearrangement.locus,
                    "species": rearrangement.locus_species,
                    "v_call": rearrangement.v_call,
                    "j_call": rearrangement.j_call,
                }

                v_domain = VDomainRef(chain_id, domain.rearrangement_index)
                binding_site = binding_sites[v_domain]
                partner = (
                    binding_site.light_domain if v_domain == binding_site.heavy_domain else binding_site.heavy_domain
                )
                row |= {
                    "binding_site_id": binding_site.binding_site_id,
                    "binding_site_type": binding_site.binding_site_type,
                    "partner_chain_id": partner.chain_id if partner else None,
                    "partner_domain_index": domain_indices[partner] if partner else None,
                    "ambiguous": binding_site.ambiguous,
                    "candidates": _format_v_domains(binding_site.candidates, domain_indices),
                }

            rows.append({field: row.get(field) for field in LAYOUT_FIELDS})
    return rows


def binding_site_rows(molecule: MoleculeAnnotation) -> list[dict[str, Any]]:
    """Binding sites of a molecule, one row per binding site."""
    domain_indices = _domain_indices(molecule)

    rows = []
    for binding_site in molecule.binding_sites:
        row: dict[str, Any] = {
            "molecule_id": molecule.molecule_id,
            "binding_site_id": binding_site.binding_site_id,
            "binding_site_type": binding_site.binding_site_type,
            "linker_length": binding_site.linker_length,
            "ambiguous": binding_site.ambiguous,
            "candidates": _format_v_domains(binding_site.candidates, domain_indices),
        }
        for prefix, v_domain in (("heavy", binding_site.heavy_domain), ("light", binding_site.light_domain)):
            if v_domain:
                rearrangement = molecule.chains[v_domain.chain_id].rearrangements[v_domain.rearrangement_index]
                row |= {
                    f"{prefix}_chain_id": v_domain.chain_id,
                    f"{prefix}_domain_index": domain_indices[v_domain],
                    f"{prefix}_locus": rearrangement.locus,
                    f"{prefix}_v_call": rearrangement.v_call,
                }
        rows.append({field: row.get(field) for field in BINDING_SITE_FIELDS})
    return rows


def _domain_indices(molecule: MoleculeAnnotation) -> dict[VDomainRef, int]:
    """1-based indices of V domains among all domains of their chains."""
    return {
        VDomainRef(chain_id, domain.rearrangement_index): domain_index
        for chain_id, chain in molecule.chains.items()
        for domain_index, domain in enumerate(chain.domains, start=1)
        if domain.domain in V_DOMAINS and domain.rearrangement_index is not None
    }


def _format_v_domains(v_domains: Sequence[VDomainRef], domain_indices: dict[VDomainRef, int]) -> Optional[str]:
    if not v_domains:
        return None
    return ";".join(f"{v_domain.chain_id}:{domain_indices[v_domain]}" for v_domain in v_domains)
