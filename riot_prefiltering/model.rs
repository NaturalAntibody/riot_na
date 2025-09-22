use pyo3::prelude::*;
use serde::Serialize;
use std::cmp::Ordering;

pub type Kmer = Vec<u8>;
pub type ByteSequence = Vec<u8>;
pub type Coverage = i32;
pub type RevComp = bool;
pub type GeneId = String;

pub struct KmerGeneIndexEntry {
    pub gene_id: GeneId,
    pub position: usize,
}
#[derive(PartialEq, Eq, Debug)]
pub struct KmerMatch {
    pub target_position: usize,
    pub query_position: usize,
}

impl PartialOrd for KmerMatch {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for KmerMatch {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.target_position == other.target_position {
            return self.query_position.cmp(&other.query_position);
        }
        self.target_position.cmp(&other.target_position)
    }
}

#[pyclass]
#[derive(PartialEq, Eq, Debug, Clone, Serialize)]
pub struct GeneSegment {
    #[pyo3(get)]
    pub start_target: usize,
    #[pyo3(get)]
    pub end_target: usize,
    #[pyo3(get)]
    pub start_query: usize,
    #[pyo3(get)]
    pub end_query: usize,
    #[pyo3(get)]
    pub coverage: Coverage,
    #[pyo3(get)]
    pub match_count: usize,
    #[pyo3(get)]
    pub gene_id: GeneId,
}

#[pymethods]
impl GeneSegment {
    fn __str__(&self) -> PyResult<String> {
        Ok(format!(
            "GeneSegment(target={}..{}, query={}..{}, coverage={}, matches={}, gene_id={})",
            self.start_target, self.end_target, self.start_query, self.end_query, self.coverage, self.match_count, self.gene_id
        ))
    }
    fn __repr__(&self) -> PyResult<String> {
        Ok(format!(
            "GeneSegment(target={}..{}, query={}..{}, coverage={}, matches={}, gene_id={})",
            self.start_target, self.end_target, self.start_query, self.end_query, self.coverage, self.match_count, self.gene_id
        ))
    }
}


#[pyclass]
#[derive(PartialEq, Eq, Debug, Clone, Serialize)]
pub struct GeneMatch {
    #[pyo3(get)]
    pub gene_id: GeneId,
    #[pyo3(get)]
    pub rev_comp: RevComp,
    #[pyo3(get)]
    pub coverage: Coverage,
}

impl PartialOrd for GeneMatch {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for GeneMatch {
    fn cmp(&self, other: &Self) -> Ordering {
        match other.coverage.cmp(&self.coverage) {
            std::cmp::Ordering::Equal => self.gene_id.cmp(&other.gene_id),
            x => x,
        }
    }
}

#[pymethods]
impl GeneMatch {
    fn __str__(&self) -> PyResult<String> {
        Ok(format!(
            "GeneMatch(gene_id={}, rev_comp={}, coverage={})",
            self.gene_id, self.rev_comp, self.coverage
        ))
    }
    fn __repr__(&self) -> PyResult<String> {
        Ok(format!(
            "GeneMatch(gene_id={}, rev_comp={}, coverage={})",
            self.gene_id, self.rev_comp, self.coverage
        ))
    }
}

#[pyclass]
pub struct PrefilteringResult {
    #[pyo3(get)]
    pub query: String,
    #[pyo3(get)]
    pub rev_comp_query: String,
    #[pyo3(get)]
    pub top_matches: Vec<GeneMatch>,
}

#[pymethods]
impl PrefilteringResult {
    fn __str__(&self) -> PyResult<String> {
        Ok(format!(
            "PrefilteringResult(query={}, rev_comp_query={}, top_matches={:?})",
            self.query, self.rev_comp_query, self.top_matches
        ))
    }
    fn __repr__(&self) -> PyResult<String> {
        Ok(format!(
            "PrefilteringResult(query={}, rev_comp_query={}, top_matches={:?})",
            self.query, self.rev_comp_query, self.top_matches
        ))
    }
}

// New segment-centric data structures

#[pyclass]
#[derive(PartialEq, Debug, Clone, Serialize)]
pub struct SegmentMatch {
    #[pyo3(get)]
    pub query_start: usize,
    #[pyo3(get)]
    pub query_end: usize,
    #[pyo3(get)]
    pub coverage: Coverage,
    #[pyo3(get)]
    pub match_count: usize,
    #[pyo3(get)]
    pub matching_genes: Vec<GeneMatch>,
    #[pyo3(get)]
    pub min_target_start: usize,
    #[pyo3(get)]
    pub max_target_start: usize,
    #[pyo3(get)]
    pub segment_start: usize,
}

impl SegmentMatch {
    pub(crate) fn recalculate_bounds_and_coverage(&mut self) {
        if self.matching_genes.is_empty() {
            self.query_start = 0;
            self.query_end = 0;
            self.coverage = 0;
            self.min_target_start = 0;
            self.max_target_start = 0;
            self.segment_start = 0;
            self.match_count = 0;
            return;
        }

        // Since GeneMatch no longer has coordinates, we only recalculate coverage
        // The segment bounds should already be set correctly when the segment is created
        let max_cov = self.matching_genes.iter().map(|gene| gene.coverage).max().unwrap_or(0);
        self.coverage = max_cov;
        self.match_count = self.matching_genes.len();
    }

}

impl PartialOrd for SegmentMatch {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for SegmentMatch {
    fn cmp(&self, other: &Self) -> Ordering {
        // Sort by query position first (for natural domain order)
        self.coverage.cmp(&other.coverage)
    }
}

impl Eq for SegmentMatch {}

#[pymethods]
impl SegmentMatch {
    fn __str__(&self) -> PyResult<String> {
        Ok(format!(
            "SegmentMatch(query={}..{}, length={}, coverage={}, matches={}, genes={})",
            self.query_start, self.query_end, self.query_length()?, self.coverage, self.match_count, self.matching_genes.len()
        ))
    }
    fn __repr__(&self) -> PyResult<String> {
        Ok(format!(
            "SegmentMatch(query={}..{}, length={}, coverage={}, matches={}, genes={:?})",
            self.query_start, self.query_end, self.query_length()?, self.coverage, self.match_count, self.matching_genes
        ))
    }

    /// Get the length of the segment in query bases
    #[pyo3(name = "query_length")]
    fn query_length(&self) -> PyResult<usize> {
        Ok(self.query_end - self.query_start)
    }
}

#[pyclass]
pub struct PrefilteringSegmentResult {
    #[pyo3(get)]
    pub query: String,
    #[pyo3(get)]
    pub rev_comp_query: String,
    #[pyo3(get)]
    pub segments: Vec<SegmentMatch>,
}

#[pymethods]
impl PrefilteringSegmentResult {
    fn __str__(&self) -> PyResult<String> {
        Ok(format!(
            "PrefilteringSegmentResult(query={}, rev_comp_query={}, segments={})",
            self.query, self.rev_comp_query, self.segments.len()
        ))
    }
    fn __repr__(&self) -> PyResult<String> {
        Ok(format!(
            "PrefilteringSegmentResult(query={}, segments={:?})",
            self.query, self.segments
        ))
    }

    /// Get the number of distinct domains detected
    #[pyo3(name = "domain_count")]
    fn domain_count(&self) -> PyResult<usize> {
        Ok(self.segments.len())
    }
}
