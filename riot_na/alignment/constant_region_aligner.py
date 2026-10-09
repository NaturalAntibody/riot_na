"""
Constant domains of whole chains.

C alignments of rearrangements cover the constant region following V (and J) only. Here constant domains are
searched for in all parts of a chain not covered by V domains, so constant domains without a preceding V domain
(e.g. Fc chains) and repeated constant domains (e.g. two Fc regions of a single chain Fc) are found as well.

Constant domains (CH1, CH2, CH3, CL) of C genes with domain definitions are aligned separately, so every domain is
assigned its own gene, which recognizes hybrid constant regions (e.g. IgG2 CH1 and hinge with IgG4 CH2 and CH3).
Hinges, short and of low complexity, are aligned only next to CH1 or CH2. C genes without domain definitions
(e.g. IGHM, IGHA, alpaca C genes) are aligned as a whole.
"""

from collections import defaultdict
from dataclasses import dataclass, replace
from typing import Optional, Sequence

from ssw_aligner import StripedSmithWaterman

from riot_na.alignment.gene_aligner import (
    AA_ALIGNER_PARAMS,
    GeneAlignerAA,
    get_aa_aligner_params,
)
from riot_na.alignment.skbio_alignment import align_aa
from riot_na.data.model import (
    AlignmentEntryAA,
    ChainType,
    GeneAA,
    GermlineGene,
    Locus,
    Organism,
)
from riot_na.data.scheme_mapping_facade import CSchemeMappingFacade
from riot_na.schemes.constant_domains import (
    HINGE,
    MEMBRANE,
    GermlineDomain,
    get_germline_domains,
)

# constant region of a C gene without domain definitions
CONSTANT_REGION = "C"
CH1 = "CH1"
CH2 = "CH2"
CH3 = "CH3"
# CHS (EU 446-447) of secreted IgG, C germlines are membrane-bound forms without it
SECRETED_TAIL = "GK"
# shortest query fragment searched for constant domains, as short as the shortest C alignment
MIN_INTERVAL_LENGTH = 4


@dataclass(frozen=True)
class QueryInterval:
    start: int
    end: int
    # C genes of this species (and of the default species) are aligned, it is the species of the preceding V domain
    species: Organism
    # C genes of this locus win ties, it is the locus of the preceding V domain
    locus: Optional[Locus] = None


@dataclass(frozen=True)
class ConstantDomainAlignment:
    domain: str
    # domain span in the query, 0-based, end exclusive. It can be longer than the alignment (by the secreted tail)
    query_start: int
    query_end: int
    # query positions refer to the whole query, germline positions to the whole germline
    alignment: AlignmentEntryAA


class ConstantRegionAligner:  # pylint: disable=too-many-instance-attributes
    def __init__(
        self,
        c_aligners: dict[Organism, dict[Locus, GeneAlignerAA]],
        imgt_c_scheme_mapping_facade: CSchemeMappingFacade,
    ) -> None:
        # constant domains not preceded by V domains are aligned to C genes of this species
        self.default_species = Organism.HOMO_SAPIENS if Organism.HOMO_SAPIENS in c_aligners else next(iter(c_aligners))

        self.genes: dict[str, GeneAA] = {}
        # species of C genes of every organism, e.g. mouse sequences are aligned to human C genes
        self.c_gene_species: dict[Organism, set[Organism]] = defaultdict(set)
        for organism, organism_aligners in c_aligners.items():
            for aligner in organism_aligners.values():
                for lookup_id, gene in aligner.gene_lookup.items():
                    assert isinstance(gene, GeneAA)
                    self.genes[lookup_id] = gene
                    self.c_gene_species[organism].add(gene.species)

        self.germline_domains: dict[str, dict[str, GermlineDomain]] = {}
        for lookup_id, gene in self.genes.items():
            imgt_scheme_mapping = imgt_c_scheme_mapping_facade.get_mapping(gene.species, gene.name)
            if imgt_scheme_mapping:
                germline_domains = get_germline_domains(
                    imgt_scheme_mapping, ChainType.from_locus(gene.locus), len(gene.sequence)
                )
                self.germline_domains[lookup_id] = {
                    germline_domain.domain: germline_domain for germline_domain in germline_domains
                }

        self.domain_aligners = self._create_domain_aligners()

    def _create_domain_aligners(self) -> dict[Organism, dict[str, GeneAlignerAA]]:
        domain_genes: dict[Organism, dict[str, list[GeneAA]]] = defaultdict(lambda: defaultdict(list))
        for lookup_id, gene in self.genes.items():
            germline_domains = self.germline_domains.get(lookup_id)
            if not germline_domains:
                domain_genes[gene.species][CONSTANT_REGION].append(gene)
                continue
            for domain, germline_domain in germline_domains.items():
                if domain not in (HINGE, MEMBRANE):
                    domain_genes[gene.species][domain].append(
                        replace(gene, sequence=gene.sequence[germline_domain.start : germline_domain.end])
                    )

        params = get_aa_aligner_params(GermlineGene.C)
        # short alignments are allowed only at the start of C genes (e.g. CH1 fragments in Fabs), so not in CH2 and CH3
        inner_domain_params = params | {
            "alignment_length_threshold": params["min_internal_alignment_length"],
            "min_internal_alignment_length": None,
        }
        return {
            species: {
                domain: GeneAlignerAA(genes, **(inner_domain_params if domain in (CH2, CH3) else params))
                for domain, genes in species_domain_genes.items()
            }
            for species, species_domain_genes in domain_genes.items()
        }

    def align_constant_domains(self, query: str, intervals: Sequence[QueryInterval]) -> list[ConstantDomainAlignment]:
        """
        Align constant domains of the query intervals. The best constant domain is taken first and the rest of its
        interval is searched again, until no constant domain is found. Then hinges are aligned next to CH1 and CH2.
        """
        domain_alignments: list[ConstantDomainAlignment] = []
        interval_alignments: dict[QueryInterval, Optional[ConstantDomainAlignment]] = {}
        pending_intervals = [interval for interval in intervals if interval.end - interval.start >= MIN_INTERVAL_LENGTH]

        while pending_intervals:
            candidates = []
            for interval in pending_intervals:
                if interval not in interval_alignments:
                    interval_alignments[interval] = self._align_interval(query, interval)
                domain_alignment = interval_alignments[interval]
                if domain_alignment is not None:
                    candidates.append((domain_alignment, interval))

            if not candidates:
                break

            domain_alignment, interval = max(candidates, key=lambda candidate: rank(candidate[0], candidate[1]))
            domain_alignments.append(domain_alignment)
            pending_intervals.remove(interval)
            pending_intervals += [
                remaining_interval
                for remaining_interval in (
                    replace(interval, end=domain_alignment.query_start),
                    replace(interval, start=domain_alignment.query_end),
                )
                if remaining_interval.end - remaining_interval.start >= MIN_INTERVAL_LENGTH
            ]

        domain_alignments += self._align_hinges(query, intervals, domain_alignments)
        return [
            self._with_secreted_tail(query, domain_alignment)
            for domain_alignment in sorted(domain_alignments, key=lambda domain_alignment: domain_alignment.query_start)
        ]

    def _align_interval(self, query: str, interval: QueryInterval) -> Optional[ConstantDomainAlignment]:
        interval_query = query[interval.start : interval.end]
        candidates = []
        for species in self.c_gene_species[interval.species] | self.c_gene_species[self.default_species]:
            for domain, aligner in self.domain_aligners[species].items():
                alignment = aligner.align(interval_query)[0].best_alignment
                if alignment is not None:
                    candidates.append(self._to_germline_alignment(query, interval.start, domain, alignment))

        return max(candidates, key=lambda candidate: rank(candidate, interval), default=None)

    def _to_germline_alignment(
        self, query: str, query_offset: int, domain: str, domain_alignment: AlignmentEntryAA
    ) -> ConstantDomainAlignment:
        gene = self.genes[domain_alignment.lookup_gene_id()]
        germline_domain = self.germline_domains.get(domain_alignment.lookup_gene_id(), {}).get(domain)
        germline_offset = germline_domain.start if germline_domain else 0

        alignment = replace(
            domain_alignment,
            q_start=domain_alignment.q_start + query_offset,
            q_end=domain_alignment.q_end + query_offset,
            t_start=domain_alignment.t_start + germline_offset,
            t_end=domain_alignment.t_end + germline_offset,
            q_seq=query,
            t_seq=gene.sequence,
        )
        return ConstantDomainAlignment(domain, alignment.q_start, alignment.q_end, alignment)

    def _align_hinges(
        self, query: str, intervals: Sequence[QueryInterval], domain_alignments: Sequence[ConstantDomainAlignment]
    ) -> list[ConstantDomainAlignment]:
        """
        Align hinges between CH1 and CH2 (or V domain and CH2, as in scFv-Fc, or CH1 and V domain, as in Fab').
        A hinge has to follow CH1 aligned up to its end or to precede CH2 aligned from its start.
        """
        hinges = []
        for interval in intervals:
            interval_domains = sorted(
                (
                    domain_alignment
                    for domain_alignment in domain_alignments
                    if interval.start <= domain_alignment.query_start < interval.end
                ),
                key=lambda domain_alignment: domain_alignment.query_start,
            )
            previous_domains: list[Optional[ConstantDomainAlignment]] = [None, *interval_domains]
            next_domains: list[Optional[ConstantDomainAlignment]] = [*interval_domains, None]
            for previous_domain, next_domain in zip(previous_domains, next_domains):
                gap_start = previous_domain.query_end if previous_domain else interval.start
                gap_end = next_domain.query_start if next_domain else interval.end
                hinge = self._align_hinge(query, gap_start, gap_end, previous_domain, next_domain)
                if hinge:
                    hinges.append(hinge)
        return hinges

    def _align_hinge(  # pylint: disable=too-many-arguments
        self,
        query: str,
        gap_start: int,
        gap_end: int,
        previous_domain: Optional[ConstantDomainAlignment],
        next_domain: Optional[ConstantDomainAlignment],
    ) -> Optional[ConstantDomainAlignment]:
        after_ch1 = previous_domain is not None and self._aligned_to_domain_end(previous_domain, CH1)
        before_ch2 = next_domain is not None and self._aligned_from_domain_start(next_domain, CH2)
        if gap_end <= gap_start or not (after_ch1 or before_ch2):
            return None

        neighbour = next_domain if before_ch2 else previous_domain
        assert neighbour is not None
        gap_query = query[gap_start:gap_end]
        aligner = StripedSmithWaterman(gap_query, **AA_ALIGNER_PARAMS)
        min_seq_identity = get_aa_aligner_params(GermlineGene.C)["min_seq_identity"]

        candidates = []
        for lookup_id, gene in self.genes.items():
            hinge = self.germline_domains.get(lookup_id, {}).get(HINGE)
            if hinge is None or (gene.species, gene.locus) != (neighbour.alignment.species, neighbour.alignment.locus):
                continue

            hinge_sequence = gene.sequence[hinge.start : hinge.end]
            alignment = align_aa(aligner, gap_query, gene.name, hinge_sequence, db_length=len(hinge_sequence))
            follows_ch1 = after_ch1 and alignment.q_start == 0 and alignment.t_start == 0
            precedes_ch2 = before_ch2 and alignment.q_end == len(gap_query) and alignment.t_end == len(hinge_sequence)
            if not (follows_ch1 or precedes_ch2) or (alignment.seq_identity or 0) < min_seq_identity:
                continue

            candidates.append(
                (
                    # the neighbour gene wins ties
                    (alignment.alignment_score or 0, gene.name == neighbour.alignment.target_id),
                    AlignmentEntryAA(
                        target_id=gene.name,
                        alignment_score=alignment.alignment_score,
                        seq_identity=alignment.seq_identity,
                        e_value=alignment.e_value,
                        q_start=gap_start + alignment.q_start,
                        q_end=gap_start + alignment.q_end,
                        t_start=hinge.start + alignment.t_start,
                        t_end=hinge.start + alignment.t_end,
                        cigar=alignment.cigar,
                        species=gene.species,
                        locus=gene.locus,
                        q_seq=query,
                        t_seq=gene.sequence,
                    ),
                )
            )

        if not candidates:
            return None
        _, hinge_alignment = max(candidates, key=lambda candidate: candidate[0])
        return ConstantDomainAlignment(HINGE, hinge_alignment.q_start, hinge_alignment.q_end, hinge_alignment)

    def _aligned_to_domain_end(self, domain_alignment: ConstantDomainAlignment, domain: str) -> bool:
        germline_domain = self.germline_domains.get(domain_alignment.alignment.lookup_gene_id(), {}).get(domain)
        return (
            domain_alignment.domain == domain
            and germline_domain is not None
            and domain_alignment.alignment.t_end == germline_domain.end
        )

    def _aligned_from_domain_start(self, domain_alignment: ConstantDomainAlignment, domain: str) -> bool:
        germline_domain = self.germline_domains.get(domain_alignment.alignment.lookup_gene_id(), {}).get(domain)
        return (
            domain_alignment.domain == domain
            and germline_domain is not None
            and domain_alignment.alignment.t_start == germline_domain.start
        )

    def _with_secreted_tail(self, query: str, domain_alignment: ConstantDomainAlignment) -> ConstantDomainAlignment:
        """
        C germlines are membrane-bound forms, so CHS of secreted IgG (GK following CH3) is not aligned.
        Add it to CH3 aligned up to its end. It is clipped later, if it overlaps the next domain.
        """
        if not self._aligned_to_domain_end(domain_alignment, CH3):
            return domain_alignment

        tail_length = 0
        if query.startswith(SECRETED_TAIL, domain_alignment.query_end):
            tail_length = len(SECRETED_TAIL)
        elif query[domain_alignment.query_end :] == SECRETED_TAIL[0]:
            # C-terminal lysine clipped
            tail_length = 1
        return replace(domain_alignment, query_end=domain_alignment.query_end + tail_length)


def rank(domain_alignment: ConstantDomainAlignment, interval: QueryInterval) -> tuple[bool, float, float, bool]:
    """
    Alignments to single domains and to whole C genes are compared, so alignments (except short C gene fragments)
    are ranked by identity before score, e.g. a human Fc is assigned human CH2 and CH3, not a whole alpaca C gene.
    """
    alignment = domain_alignment.alignment
    min_internal_alignment_length = get_aa_aligner_params(GermlineGene.C)["min_internal_alignment_length"]
    return (
        alignment.t_end - alignment.t_start >= min_internal_alignment_length,
        alignment.seq_identity or 0,
        alignment.alignment_score or 0,
        alignment.locus == interval.locus,
    )
