class GeneMatch:
    gene_id: str
    rev_comp: bool
    coverage: int

class GeneSegment:
    start_target: int
    end_target: int
    start_query: int
    end_query: int
    coverage: int
    match_count: int

class SegmentMatch:
    query_start: int
    query_end: int
    coverage: int
    match_count: int
    matching_genes: list[GeneMatch]
    segment_start: int

    def query_length(self) -> int: ...

class PrefilteringResult:
    query: str
    rev_comp_query: str
    top_matches: list[GeneMatch]

class PrefilteringSegmentResult:
    query: str
    rev_comp_query: str
    segments: list[SegmentMatch]

    def domain_count(self) -> int: ...

class Prefiltering:
    def __init__(
        self,
        genes: dict[str, str],
        kmer_size: int,
        distance_threshold: int,
        top_n: int,
        modulo_n: int,
        min_segment_length: int = 30,
        min_coverage: int = 20,
    ): ...
    def calculate_top_matches_with_rev_comp(self, query: str) -> PrefilteringResult: ...
    def calculate_top_matches(self, query: str) -> PrefilteringResult: ...
    def calculate_segment_matches(self, query: str) -> PrefilteringSegmentResult: ...
    def calculate_segment_matches_with_rev_comp(self, query: str) -> PrefilteringSegmentResult: ...
