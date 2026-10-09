"""
Domains of C germlines (CH1, hinge, CH2, CH3, CL), derived from their IMGT C scheme mappings: IMGT positions
of C germlines are prefixed with their domain (CH1_1.4, H_1, H1_1, CH2_84.4, CH3_125, CL_1.4).
"""

from dataclasses import dataclass
from itertools import groupby

from riot_na.data.model import AlignmentString, ChainType, CScheme
from riot_na.data.scheme_definitions import get_c_scheme_positions

HINGE = "H"
# germline residues past the numbered domains (membrane exons)
MEMBRANE = "M"


@dataclass(frozen=True)
class GermlineDomain:
    domain: str
    # 0-based germline positions, end exclusive
    start: int
    end: int


def domain_of_imgt_position(imgt_position: str) -> str:
    domain = imgt_position.split("_")[0]
    # IGHG3 hinge is encoded by 4 exons: H1-H4
    return HINGE if domain.startswith(HINGE) else domain


def get_germline_domains(
    imgt_scheme_mapping: AlignmentString, chain_type: ChainType, germline_length: int
) -> list[GermlineDomain]:
    imgt_positions = get_c_scheme_positions(CScheme.IMGT, chain_type)

    domains: list[str] = []
    position_index = 0
    for op in imgt_scheme_mapping:
        if op == "M":
            domains.append(domain_of_imgt_position(imgt_positions[position_index]))
            position_index += 1
        elif op == "D":
            position_index += 1
        elif op == "I":
            # IMGT insertions are inside domains, except 1.x which precede position 1 of the domain,
            # so an insertion belongs to the domain of the next scheme position
            domains.append(domain_of_imgt_position(imgt_positions[min(position_index, len(imgt_positions) - 1)]))
    domains += [MEMBRANE] * (germline_length - len(domains))

    germline_domains = []
    for domain, positions in groupby(enumerate(domains), key=lambda position_domain: position_domain[1]):
        domain_positions = [position for position, _ in positions]
        germline_domains.append(GermlineDomain(domain, domain_positions[0], domain_positions[-1] + 1))
    return germline_domains
