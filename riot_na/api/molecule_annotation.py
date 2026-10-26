"""
Molecule annotation: binding sites formed by V domains of the chains of a molecule.

V domains are paired by the format of the chains, in this order:
- scFv: adjacent heavy and light V domains of a chain (no constant domain in between) joined by a linker of at least
  MIN_SCFV_LINKER_LENGTH residues, unless one of them is given in pairs as a V domain with other V domains,
- diabody: the same with a shorter linker, too short to pair V domains of one chain. V domains pair with V domains of
  another copy of the chain (homodimer) or, if two chains are each other's only candidates, of the other chain
  (heterodimer). Units of one chain pair with each other only if they are given in pairs with each other (e.g. tandem
  diabody),
- Fab: V domain followed by CH1 pairs with V domain followed by CL, if one of them is heavy and the other light
  (which pairs domain crossovers of CrossMabs too). Constant regions of C genes without domain definitions
  (e.g. IGHM) aligned from their start are CH1 of heavy and CL of light C genes. V domains stacked in front of them
  (e.g. outer V domains of DVD-Ig) pair by their position. Fab arms of one chain pair only if joined by a linker of at
  least MIN_SCFAB_LINKER_LENGTH residues (single chain Fab), then only with each other. Fab arms with more than one
  possible partner are ambiguous,
- Fv: V domains given in pairs (pairing hints), which are not paired by their format (e.g. V domain only chains).
  A light V domain given in pairs with many heavy V domains pairs with each of them (common light chain),
- single domain: other heavy V domains (e.g. VHH), other light V domains and V domains of Fab arms are unpaired.

Pairing hints (pairs) are pairs of chains or V domains of chains which pair with each other. V domains given in pairs
pair only with V domains they are given with, which resolves ambiguous Fab arms and diabodies. Light chain Fab arms
(followed by CL) and light V domains paired as Fv are not reserved by pairs, they can be shared by heavy chain Fab arms
or heavy V domains (common light chain). Pairs choose partners of diabody units. Chains given in pairs (chain ID
only) don't prevent scFv of their V domains, only V domains given in pairs (chain ID and domain index) do.
"""

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Mapping, Optional, Sequence, Union

from riot_na.alignment.constant_region_aligner import CONSTANT_REGION
from riot_na.data.model import (
    BindingSite,
    ChainAnnotation,
    ChainDomain,
    Locus,
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
# constant regions of C genes without domain definitions aligned from (about) the C gene start contain CH1 or CL
MAX_CONSTANT_REGION_GERMLINE_START = 3

FAB = "Fab"
FV = "Fv"
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

# pairing hint: chain ID (all V domains of the chain) or chain ID and domain index of a V domain
# (1-based position among all domains of the chain, as in molecule layout and binding site rows)
PairingTarget = Union[str, tuple[str, int]]
# chain ID and position of a V domain in the chain domains, None for all V domains of the chain
_Target = tuple[str, Optional[int]]


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

    def matches(self, target: _Target) -> bool:
        chain_id, position = target
        return self.chain_id == chain_id and position in (None, self.position)


@dataclass(frozen=True)
class _FabArm:
    # V domains stacked in front of the constant domain, the last one is followed by it
    v_domains: tuple[_VDomain, ...]
    constant_domain: ChainDomain = field(compare=False)
    # CH1 or CL
    constant_domain_type: str = field(compare=False)

    @property
    def chain_id(self) -> str:
        return self.v_domains[-1].chain_id

    @property
    def heavy(self) -> bool:
        return self.v_domains[-1].heavy


class _BindingSitePairing:  # pylint: disable=too-few-public-methods
    def __init__(
        self,
        chains: Mapping[str, ChainAnnotation],
        residue_length: int,
        pairs: Optional[Sequence[tuple[PairingTarget, PairingTarget]]],
    ) -> None:
        self.chains = chains
        self.residue_length = residue_length
        self.v_domains = {
            chain_id: [
                _VDomain(chain_order, chain_id, position, domain)
                for position, domain in enumerate(chain.domains)
                if domain.domain in V_DOMAINS
            ]
            for chain_order, (chain_id, chain) in enumerate(chains.items())
        }
        self.hints = [(self._target(first), self._target(second)) for first, second in pairs or ()]
        self.binding_sites: list[BindingSite] = []
        self.paired: set[_VDomain] = set()
        # V domains of Fab arms are unpaired, if no partner is found
        self.fab_arm_v_domains: set[_VDomain] = set()

    def _target(self, target: PairingTarget) -> _Target:
        chain_id, domain_index = (target, None) if isinstance(target, str) else target
        if chain_id not in self.chains:
            raise ValueError(f"Unknown chain in pairs: {chain_id}")
        if domain_index is None:
            return chain_id, None
        if all(v_domain.position != domain_index - 1 for v_domain in self.v_domains[chain_id]):
            raise ValueError(f"No V domain with domain index {domain_index} in chain {chain_id} in pairs")
        return chain_id, domain_index - 1

    def _hint_partners(self, v_domains: Sequence[_VDomain]) -> list[_Target]:
        """Targets which any of the V domains is given in pairs with."""
        return [
            partner
            for hint in self.hints
            for target, partner in (hint, hint[::-1])
            if any(v_domain.matches(target) for v_domain in v_domains)
        ]

    def _hinted(self, v_domains: Sequence[_VDomain], partner_v_domains: Sequence[_VDomain]) -> bool:
        """V domains are given in pairs with partner V domains."""
        return any(
            partner_v_domain.matches(target)
            for target in self._hint_partners(v_domains)
            for partner_v_domain in partner_v_domains
        )

    def _hints_allow(self, v_domains: Sequence[_VDomain], partner_v_domains: Sequence[_VDomain]) -> bool:
        """V domains given in pairs pair only with V domains they are given with."""
        hinted = (self._hint_partners(v_domains), self._hint_partners(partner_v_domains))
        return not any(hinted) or self._hinted(v_domains, partner_v_domains)

    def _domain_hints_allow(self, first: _VDomain, second: _VDomain) -> bool:
        """
        V domains given in pairs as V domains (chain ID and domain index) pair only with V domains they are given with.
        """
        return all(
            not partners or any(other.matches(partner) for partner in partners)
            for v_domain, other in ((first, second), (second, first))
            for partners in [
                [
                    partner
                    for hint in self.hints
                    for target, partner in (hint, hint[::-1])
                    if target == (v_domain.chain_id, v_domain.position)
                ]
            ]
        )

    def pair(self) -> list[BindingSite]:
        self._pair_single_chain_fvs()
        self._pair_fab_arms()
        self._pair_hinted_v_domains()

        for chain_v_domains in self.v_domains.values():
            for v_domain in chain_v_domains:
                if v_domain not in self.paired:
                    single_domain = v_domain.heavy and v_domain not in self.fab_arm_v_domains
                    self._add_binding_site(SINGLE_DOMAIN if single_domain else UNPAIRED, [v_domain])

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
        diabody_units: list[tuple[_VDomain, _VDomain]] = []

        for chain_v_domains in self.v_domains.values():
            i = 0
            while i + 1 < len(chain_v_domains):
                first, second = chain_v_domains[i], chain_v_domains[i + 1]
                if second.position != first.position + 1 or first.heavy == second.heavy:
                    i += 1
                    continue

                linker_length = self._linker_length(first.domain, second.domain)
                if linker_length < MIN_SCFV_LINKER_LENGTH:
                    diabody_units.append((first, second))
                    self.paired.update((first, second))
                elif self._domain_hints_allow(first, second):
                    self._add_binding_site(SCFV, [first, second], linker_length=linker_length)
                else:
                    # V domains given in pairs (as V domains) with other V domains don't form a scFv
                    i += 1
                    continue
                i += 2

        self._pair_diabodies(diabody_units)

    def _pair_diabodies(self, units: Sequence[tuple[_VDomain, _VDomain]]) -> None:
        """
        Diabody units (V domains joined by short linkers) of chains which are each other's only candidates
        form a heterodimer, other units pair with another copy of their chain (homodimer). Units of one chain are
        candidates of each other only if they are given in pairs with each other (e.g. tandem diabody). V domains of
        units which can't pair with any unit, because their V domains are given in pairs (as V domains) with other
        V domains, are paired by pairs.
        """
        candidates = {
            unit: [
                other_unit
                for other_unit in units
                if (
                    other_unit[0].chain_id != unit[0].chain_id
                    or (other_unit != unit and self._hinted(unit, other_unit))
                )
                and self._hints_allow(unit, other_unit)
            ]
            for unit in units
        }
        heterodimer_units: set[tuple[_VDomain, _VDomain]] = set()

        for unit in units:
            first, second = unit
            unit_candidates = candidates[unit]
            if unit in heterodimer_units:
                continue
            if len(unit_candidates) == 1 and candidates[unit_candidates[0]] == [unit]:
                other_first, other_second = unit_candidates[0]
                heterodimer_units.update((unit, unit_candidates[0]))
                # every V domain pairs with the complementary V domain of the other chain
                self._add_binding_site(
                    DIABODY, [first, other_second if first.heavy != other_second.heavy else other_first]
                )
                self._add_binding_site(
                    DIABODY, [second, other_first if second.heavy != other_first.heavy else other_second]
                )
            elif not unit_candidates and not self._domain_hints_allow(first, second):
                self.paired.difference_update(unit)
            else:
                self._add_binding_site(
                    DIABODY,
                    [first, second],
                    ambiguous=bool(unit_candidates),
                    candidates=[v_domain for other_unit in unit_candidates for v_domain in other_unit],
                    linker_length=self._linker_length(first.domain, second.domain),
                )

    def _fab_arms(self) -> list[_FabArm]:
        """V domains followed by CH1 or CL, with V domains stacked in front of them."""
        fab_arms = []
        for chain_id, chain in self.chains.items():
            domains = chain.domains
            v_domains = {v_domain.position: v_domain for v_domain in self.v_domains[chain_id]}

            for position, domain in enumerate(domains):
                v_domain = v_domains.get(position - 1)
                constant_domain_type = _fab_constant_domain_type(chain, domain)
                full_constant_domain = (
                    self._domain_length(domain) >= MIN_FAB_CONSTANT_DOMAIN_LENGTH or position == len(domains) - 1
                )
                if constant_domain_type is None or not full_constant_domain:
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

                fab_arms.append(_FabArm(tuple(stack), domain, constant_domain_type))
        return fab_arms

    def _fab_partners(self, fab_arm: _FabArm, fab_arms: Sequence[_FabArm]) -> list[_FabArm]:
        partners = []
        for other_arm in fab_arms:
            if (
                other_arm.constant_domain_type == fab_arm.constant_domain_type
                or other_arm.heavy == fab_arm.heavy
                or len(other_arm.v_domains) != len(fab_arm.v_domains)
            ):
                continue
            # only CH1 arms are restricted by pairs, CL arms can be shared by many CH1 arms (common light chain)
            if self._hint_partners(fab_arm.v_domains) and not self._hinted(fab_arm.v_domains, other_arm.v_domains):
                continue
            if other_arm.chain_id == fab_arm.chain_id and not self._hinted(fab_arm.v_domains, other_arm.v_domains):
                first, second = sorted((fab_arm, other_arm), key=lambda arm: arm.constant_domain.sequence_start)
                # single chain Fab
                if self._linker_length(first.constant_domain, second.v_domains[0].domain) < MIN_SCFAB_LINKER_LENGTH:
                    continue
            partners.append(other_arm)
        return partners

    def _pair_fab_arms(self) -> None:
        fab_arms = self._fab_arms()
        for fab_arm in fab_arms:
            self.fab_arm_v_domains.update(fab_arm.v_domains)
        ch1_arms = [fab_arm for fab_arm in fab_arms if fab_arm.constant_domain_type == CH1]
        cl_arms = [fab_arm for fab_arm in fab_arms if fab_arm.constant_domain_type == CL]
        fab_arm_partners = {fab_arm: self._fab_partners(fab_arm, cl_arms) for fab_arm in ch1_arms}
        # Fab arms of a single chain Fab are tethered by its linker, they pair only with each other
        single_chain_partners = {
            fab_arm: [partner for partner in partners if partner.chain_id == fab_arm.chain_id]
            for fab_arm, partners in fab_arm_partners.items()
        }
        tethered_arms = {partner for partners in single_chain_partners.values() for partner in partners}
        fab_arm_partners = {
            fab_arm: single_chain_partners[fab_arm] or [partner for partner in partners if partner not in tethered_arms]
            for fab_arm, partners in fab_arm_partners.items()
        }
        fab_arm_partners = _one_to_one_partners(ch1_arms, cl_arms, fab_arm_partners) or fab_arm_partners
        cl_arms_with_partners: set[_FabArm] = set()

        for fab_arm in ch1_arms:
            partners = fab_arm_partners[fab_arm]
            cl_arms_with_partners.update(partners)
            stack_size = len(fab_arm.v_domains)

            for position, v_domain in enumerate(fab_arm.v_domains):
                binding_site_type = FAB if position == stack_size - 1 else FV
                if len(partners) == 1:
                    self._add_binding_site(binding_site_type, [v_domain, partners[0].v_domains[position]])
                elif partners:
                    self._add_binding_site(
                        binding_site_type,
                        [v_domain],
                        ambiguous=True,
                        candidates=[partner.v_domains[position] for partner in partners],
                    )

        # V domains of CL arms with ambiguous partners are listed as candidates of the CH1 arms
        for fab_arm in cl_arms_with_partners:
            self.paired.update(fab_arm.v_domains)

    def _pair_hinted_v_domains(self) -> None:
        """
        Pair V domains given in pairs, which are not paired by their format (e.g. V domain only chains). Light V domains
        are not reserved, as light chain Fab arms: a light V domain given in pairs with many heavy V domains pairs with
        each of them (common light chain).
        """
        unpaired_v_domains = [
            v_domain
            for chain_v_domains in self.v_domains.values()
            for v_domain in chain_v_domains
            if v_domain not in self.paired
        ]

        for v_domain in unpaired_v_domains:
            if not v_domain.heavy:
                continue
            partners = [
                partner for partner in unpaired_v_domains if not partner.heavy and self._hinted([v_domain], [partner])
            ]
            if len(partners) == 1:
                self._add_binding_site(FV, [v_domain, partners[0]])
            elif partners:
                self._add_binding_site(FV, [v_domain], ambiguous=True, candidates=partners)
                self.paired.update(partners)


def _fab_constant_domain_type(chain: ChainAnnotation, domain: ChainDomain) -> Optional[str]:
    """CH1 or CL, if the domain is a Fab constant domain."""
    if domain.domain in FAB_CONSTANT_DOMAINS:
        return domain.domain
    if domain.domain != CONSTANT_REGION or domain.constant_domain_index is None:
        return None
    constant_domain = chain.constant_domains[domain.constant_domain_index]
    if constant_domain.germline_start_aa > MAX_CONSTANT_REGION_GERMLINE_START:
        return None
    return CH1 if constant_domain.locus == Locus.IGH.value else CL


def _one_to_one_partners(
    ch1_arms: Sequence[_FabArm], cl_arms: Sequence[_FabArm], fab_arm_partners: Mapping[_FabArm, list[_FabArm]]
) -> Optional[dict[_FabArm, list[_FabArm]]]:
    """
    If there are as many CH1 as CL Fab arms, every CL arm likely pairs with one CH1 arm (e.g. FIT-Ig). Partners
    of CH1 arms are then the CL arms they are paired with in any one-to-one pairing. Returns None, if the numbers
    of arms differ (e.g. common light chain pairs with many heavy chains) or there is no one-to-one pairing.

    One-to-one pairings are perfect matchings of the arms, so a CL arm is paired with a CH1 arm in some one-to-one
    pairing, if the CH1 arm is reachable from the CL arm by an alternating path: from a CL arm to its CH1 arm in a
    found one-to-one pairing, then to another partner of this CH1 arm and so on. Arms of the path closed by the CH1
    arm can exchange their partners, which gives another one-to-one pairing.
    """
    if len(ch1_arms) != len(cl_arms):
        return None

    pairing = _one_to_one_pairing(ch1_arms, fab_arm_partners)
    if pairing is None:
        return None
    paired_ch1_arms = {cl_arm: ch1_arm for ch1_arm, cl_arm in pairing.items()}

    def reachable_ch1_arms(cl_arm: _FabArm) -> set[_FabArm]:
        reached: set[_FabArm] = set()
        pending = [paired_ch1_arms[cl_arm]]
        while pending:
            ch1_arm = pending.pop()
            if ch1_arm not in reached:
                reached.add(ch1_arm)
                pending += [paired_ch1_arms[partner] for partner in fab_arm_partners[ch1_arm]]
        return reached

    reachable = {cl_arm: reachable_ch1_arms(cl_arm) for cl_arm in cl_arms}
    return {
        fab_arm: [partner for partner in fab_arm_partners[fab_arm] if fab_arm in reachable[partner]]
        for fab_arm in ch1_arms
    }


def _one_to_one_pairing(
    ch1_arms: Sequence[_FabArm], fab_arm_partners: Mapping[_FabArm, list[_FabArm]]
) -> Optional[dict[_FabArm, _FabArm]]:
    """A one-to-one pairing of all CH1 arms with their partners (maximum bipartite matching), None if there is none."""
    paired_ch1_arms: dict[_FabArm, _FabArm] = {}

    def pair(ch1_arm: _FabArm, visited: set[_FabArm]) -> bool:
        # pair the CH1 arm with a free partner or with a partner whose CH1 arm can be paired with another partner
        for cl_arm in fab_arm_partners[ch1_arm]:
            if cl_arm not in visited:
                visited.add(cl_arm)
                if cl_arm not in paired_ch1_arms or pair(paired_ch1_arms[cl_arm], visited):
                    paired_ch1_arms[cl_arm] = ch1_arm
                    return True
        return False

    if not all(pair(ch1_arm, set()) for ch1_arm in ch1_arms):
        return None
    return {ch1_arm: cl_arm for cl_arm, ch1_arm in paired_ch1_arms.items()}


def pair_binding_sites(
    chains: Mapping[str, ChainAnnotation],
    residue_length: int = 1,
    pairs: Optional[Sequence[tuple[PairingTarget, PairingTarget]]] = None,
) -> list[BindingSite]:
    """
    Binding sites of the V domains of a molecule chains.

    - residue_length: length of residues in chain sequence positions, 3 for nucleotide chains.
    - pairs: pairing hints, pairs of chains (chain ID) or V domains (chain ID and domain index) which pair.
    V domains given in pairs pair only with V domains they are given with, except light chain Fab arms (followed by
    CL) and light V domains paired as Fv, which can be shared by many heavy chain Fab arms or heavy V domains (common
    light chain). They resolve ambiguous Fab arms (e.g. of bispecific antibodies with two different heavy and light
    chains, given as chains, or tandem Fabs with two different light chains, given as V domains) and diabodies, and
    pair V domains not paired by their format (e.g. V domain only chains) as Fv.
    """
    return _BindingSitePairing(chains, residue_length, pairs).pair()


def molecule_layout_rows(molecule: MoleculeAnnotation) -> list[dict[str, Any]]:
    """
    Domains of all chains of a molecule, one row per domain, with binding sites of V domains. V domains of more than
    one binding site (e.g. common light chains) have values of all of them, separated by ";". V domains which are
    only candidates of ambiguous binding sites (e.g. light chains of bispecific antibodies without pairing hints) have
    these binding sites and their V domains as candidates.
    """
    domain_indices = _domain_indices(molecule)
    member_binding_sites: dict[VDomainRef, list[BindingSite]] = defaultdict(list)
    candidate_binding_sites: dict[VDomainRef, list[BindingSite]] = defaultdict(list)
    for binding_site in molecule.binding_sites:
        for v_domain in _binding_site_v_domains(binding_site):
            member_binding_sites[v_domain].append(binding_site)
        for candidate in binding_site.candidates:
            candidate_binding_sites[candidate].append(binding_site)

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
                row |= _v_domain_binding_sites(
                    v_domain, member_binding_sites[v_domain], candidate_binding_sites[v_domain], domain_indices
                )

            rows.append({field: row.get(field) for field in LAYOUT_FIELDS})
    return rows


def _v_domain_binding_sites(
    v_domain: VDomainRef,
    member_binding_sites: Sequence[BindingSite],
    candidate_binding_sites: Sequence[BindingSite],
    domain_indices: dict[VDomainRef, int],
) -> dict[str, Any]:
    """
    Layout columns of binding sites of a V domain: binding sites the V domain is in, with their partner V domains
    and candidates, or else ambiguous binding sites the V domain is a candidate of, with their V domains as candidates.
    """
    if member_binding_sites:
        binding_sites = member_binding_sites
        partners = [
            partner
            for binding_site in binding_sites
            for partner in _binding_site_v_domains(binding_site)
            if partner != v_domain
        ]
        candidates = [candidate for binding_site in binding_sites for candidate in binding_site.candidates]
    else:
        binding_sites = candidate_binding_sites
        partners = []
        candidates = [
            candidate for binding_site in binding_sites for candidate in _binding_site_v_domains(binding_site)
        ]

    return {
        "binding_site_id": _join_values([binding_site.binding_site_id for binding_site in binding_sites]),
        "binding_site_type": _join_values([binding_site.binding_site_type for binding_site in binding_sites]),
        "partner_chain_id": _join_values([partner.chain_id for partner in partners]),
        "partner_domain_index": _join_values([domain_indices[partner] for partner in partners]),
        # binding sites with candidates are ambiguous
        "ambiguous": any(binding_site.ambiguous for binding_site in binding_sites),
        "candidates": _format_v_domains(list(dict.fromkeys(candidates)), domain_indices),
    }


def _binding_site_v_domains(binding_site: BindingSite) -> list[VDomainRef]:
    return [v_domain for v_domain in (binding_site.heavy_domain, binding_site.light_domain) if v_domain]


def _join_values(values: Sequence[Any]) -> Any:
    """The value, values separated by ";" if there are more of them, None if there are none."""
    if len(values) > 1:
        return ";".join(str(value) for value in values)
    return values[0] if values else None


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
