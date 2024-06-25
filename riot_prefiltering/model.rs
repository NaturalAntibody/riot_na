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
