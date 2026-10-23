"""
Chain annotation: domains of a chain in sequence order, rearrangements (V domains) and constant domains.
"""

from bisect import bisect_left
from dataclasses import dataclass, replace
from typing import Optional, Sequence

from riot_na.alignment.constant_region_aligner import (
    ConstantDomainAlignment,
    ConstantRegionAligner,
    QueryInterval,
)
from riot_na.alignment.translate import translate
from riot_na.data.model import ChainDomain, ConstantDomain, CScheme, Locus, Organism
from riot_na.schemes.scheme_alignment import SchemeAligner

V_DOMAIN_NAMES = {Locus.IGH: "VH", Locus.IGK: "VK", Locus.IGL: "VL"}
FRAMES = (0, 1, 2)
SPLIT_CODON_OVERLAP = 2
STOP_CODON = "*"
COMPLEMENT = str.maketrans("ACGTURYKMBVDHSWN", "TGCAAYRMKVBHDSWN")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


@dataclass(frozen=True)
class VDomainSpan:
    """V domain (V and J) of a chain, 0-based positions in the chain sequence, end exclusive."""

    start: int
    end: int
    species: Organism
    locus: Locus


@dataclass(frozen=True)
class TranslatedDomainAlignment:
    """Constant domain aligned to the translation of a chain (amino acid chains are their own translations)."""

    translation: str
    # reading frame of nucleotide chains, None for amino acid chains
    frame: Optional[int]
    domain_alignment: ConstantDomainAlignment

    @property
    def start(self) -> int:
        return self.chain_position(self.domain_alignment.query_start)

    @property
    def end(self) -> int:
        return self.chain_position(self.domain_alignment.query_end)

    def chain_position(self, position: int) -> int:
        return position if self.frame is None else self.frame + 3 * position

    def translation_position(self, chain_position: int) -> int:
        return chain_position if self.frame is None else (chain_position - self.frame) // 3


def constant_region_intervals(
    sequence_length: int, v_domains: Sequence[VDomainSpan], default_species: Organism
) -> list[QueryInterval]:
    """
    Parts of a chain not covered by V domains. C genes of the species of the preceding V domain are aligned to them,
    C genes of the default species if no V domain precedes them.
    """
    intervals = []
    start, species, locus = 0, default_species, None
    for v_domain in v_domains:
        if v_domain.start > start:
            intervals.append(QueryInterval(start, v_domain.start, species, locus, follows_v_domain=locus is not None))
        start, species, locus = max(start, v_domain.end), v_domain.species, v_domain.locus

    if sequence_length > start:
        intervals.append(QueryInterval(start, sequence_length, species, locus, follows_v_domain=locus is not None))
    return intervals


def split_at_stop_codons(intervals: Sequence[QueryInterval], translation: str) -> list[QueryInterval]:
    """
    Constant domains can't span stop codons, e.g. in wrong reading frames of nucleotide chains.
    """
    split_intervals = []
    for interval in intervals:
        start = interval.start
        for position in range(interval.start, interval.end):
            if translation[position] == STOP_CODON:
                split_intervals.append(replace(interval, start=start, end=position))
                start = position + 1
                interval = replace(interval, follows_v_domain=False)
        split_intervals.append(replace(interval, start=start))
    return [interval for interval in split_intervals if interval.end > interval.start]


def align_constant_domains_aa(
    constant_region_aligner: ConstantRegionAligner, sequence: str, v_domains: Sequence[VDomainSpan]
) -> list[TranslatedDomainAlignment]:
    intervals = constant_region_intervals(len(sequence), v_domains, constant_region_aligner.default_species)
    return [
        TranslatedDomainAlignment(sequence, None, domain_alignment)
        for domain_alignment in constant_region_aligner.align_constant_domains(
            sequence, split_at_stop_codons(intervals, sequence)
        )
    ]


def align_constant_domains_nt(
    constant_region_aligner: ConstantRegionAligner, sequence: str, v_domains: Sequence[VDomainSpan]
) -> list[TranslatedDomainAlignment]:
    """
    Align constant domains of a nucleotide chain in all reading frames. Constant domains found in different frames
    can overlap, the better scoring ones are kept.
    """
    # J germlines end with the first nucleotide of the first C codon (split by the J-C splice site),
    # so constant domains can start up to 2 nucleotides before the end of V domains
    intervals = constant_region_intervals(
        len(sequence),
        [replace(v_domain, end=max(v_domain.start, v_domain.end - SPLIT_CODON_OVERLAP)) for v_domain in v_domains],
        constant_region_aligner.default_species,
    )

    frame_domain_alignments = []
    for frame in FRAMES:
        translation = translate(sequence, frame)
        # residues encoded entirely within the intervals
        frame_intervals = split_at_stop_codons(
            [
                replace(interval, start=-(-(interval.start - frame) // 3), end=(interval.end - frame) // 3)
                for interval in intervals
            ],
            translation,
        )
        frame_domain_alignments += [
            TranslatedDomainAlignment(translation, frame, domain_alignment)
            for domain_alignment in constant_region_aligner.align_constant_domains(translation, frame_intervals)
        ]

    selected: list[TranslatedDomainAlignment] = []
    for domain_alignment in sorted(
        frame_domain_alignments,
        key=lambda frame_domain_alignment: frame_domain_alignment.domain_alignment.alignment.alignment_score or 0,
        reverse=True,
    ):
        if all(domain_alignment.end <= other.start or domain_alignment.start >= other.end for other in selected):
            selected.append(domain_alignment)
    return sorted(selected, key=lambda domain_alignment: domain_alignment.start)


def chain_domains(v_domains: Sequence[VDomainSpan], constant_domains: Sequence[ConstantDomain]) -> list[ChainDomain]:
    """V and constant domains of a chain in sequence order."""
    domains = [
        ChainDomain(V_DOMAIN_NAMES[v_domain.locus], v_domain.start + 1, v_domain.end, rearrangement_index=i)
        for i, v_domain in enumerate(v_domains)
    ] + [
        ChainDomain(
            constant_domain.domain,
            constant_domain.sequence_start,
            constant_domain.sequence_end,
            rearrangement_index=constant_domain.rearrangement_index,
            constant_domain_index=i,
        )
        for i, constant_domain in enumerate(constant_domains)
    ]
    return sorted(domains, key=lambda domain: (domain.sequence_start, domain.sequence_end))


def chain_layout(sequence_length: int, domains: Sequence[ChainDomain]) -> str:
    """
    Domains of a chain in sequence order, e.g. VH-CH1-H-CH2-CH3-[15]-VH-[20]-VK, [n] - n unannotated residues
    (nucleotides for nucleotide chains).
    """
    layout = []
    position = 0
    for domain in domains:
        if domain.sequence_start - 1 > position:
            layout.append(f"[{domain.sequence_start - 1 - position}]")
        layout.append(domain.domain)
        position = max(position, domain.sequence_end)

    if sequence_length > position:
        layout.append(f"[{sequence_length - position}]")
    return "-".join(layout)


def annotate_constant_domains(
    scheme_aligner: SchemeAligner,
    domain_alignments: Sequence[TranslatedDomainAlignment],
    v_domains: Sequence[VDomainSpan],
    c_scheme: Optional[CScheme],
) -> list[ConstantDomain]:
    """
    Constant domains of domain alignments, numbered in c_scheme. Every constant domain is assigned
    the rearrangement (V domain) it follows.
    """
    v_domain_starts = [v_domain.start for v_domain in v_domains]
    annotated_domains = []

    for i, translated_domain_alignment in enumerate(domain_alignments):
        domain_alignment = translated_domain_alignment.domain_alignment
        alignment = domain_alignment.alignment
        translation = translated_domain_alignment.translation

        sequence_start = translated_domain_alignment.chain_position(domain_alignment.query_start)
        sequence_end = translated_domain_alignment.chain_position(domain_alignment.query_end)
        # the secreted tail added to CH3 can't overlap the next domain
        next_v_domain = bisect_left(v_domain_starts, sequence_start)
        if next_v_domain < len(v_domains):
            sequence_end = min(sequence_end, v_domains[next_v_domain].start)
        if i + 1 < len(domain_alignments):
            next_domain = domain_alignments[i + 1]
            sequence_end = min(sequence_end, next_domain.chain_position(next_domain.domain_alignment.query_start))
        query_end = translated_domain_alignment.translation_position(sequence_end)

        constant_domain = ConstantDomain(
            domain=domain_alignment.domain,
            c_call=alignment.target_id,
            locus=alignment.locus.value,
            species=alignment.species.value,
            sequence_start=sequence_start + 1,
            sequence_end=sequence_end,
            sequence_aa=translation[domain_alignment.query_start : query_end],
            germline_start_aa=alignment.t_start + 1,
            germline_end_aa=alignment.t_end,
            cigar_aa=alignment.cigar,
            identity_aa=alignment.seq_identity or 0.0,
            rearrangement_index=next_v_domain - 1 if next_v_domain > 0 else None,
        )

        if c_scheme:
            c_numbering = scheme_aligner.number_c_alignment(translation, alignment, c_scheme)
            if c_numbering:
                c_scheme_residue_mapping, c_positional_scheme_mapping = c_numbering
                constant_domain.c_numbering_scheme = c_scheme.value
                constant_domain.c_scheme_residue_mapping = c_scheme_residue_mapping
                constant_domain.c_positional_scheme_mapping = {
                    translated_domain_alignment.chain_position(position): scheme_position
                    for position, scheme_position in c_positional_scheme_mapping.items()
                }

        annotated_domains.append(constant_domain)

    return annotated_domains
