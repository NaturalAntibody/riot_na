class GeneMatch:
    gene_id: str
    rev_comp: bool
    coverage: int

class PrefilteringResult:
    query: str
    rev_comp_query: str
    top_matches: list[GeneMatch]

class Prefiltering:
    def __init__(self, genes: dict[str, str], kmer_size: int, distance_threshold: int, top_n: int, modulo_n: int): ...
    def calculate_top_matches_with_rev_comp(self, query: str) -> PrefilteringResult: ...
    def calculate_top_matches(self, query: str) -> PrefilteringResult: ...
