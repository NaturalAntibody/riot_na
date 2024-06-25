use pyo3::{exceptions::PyValueError, prelude::*};

use crate::model::{
    ByteSequence, Coverage, GeneId, GeneMatch, Kmer, KmerGeneIndexEntry, KmerMatch,
    PrefilteringResult,
};

use ahash::{AHashMap, RandomState};
use bio::alphabets::dna;
use itertools::Itertools;
use std::collections::HashMap;

const RANDOM_SEED: usize = 42;

#[pyclass(frozen, subclass)]
pub struct Prefiltering {
    kmer_size: usize,
    distance_threshold: i32,
    top_n: usize,
    modulo_n: usize,
    kmer_genes_lookup: AHashMap<Kmer, Vec<KmerGeneIndexEntry>>,
}

#[pymethods]
impl Prefiltering {
    #[new]
    pub fn new(
        genes: HashMap<String, String>,
        kmer_size: usize,
        distance_threshold: i32,
        top_n: usize,
        modulo_n: usize,
    ) -> Self {
        let mut prefiltering: Prefiltering = Self {
            kmer_size,
            distance_threshold,
            kmer_genes_lookup: AHashMap::with_hasher(RandomState::with_seed(RANDOM_SEED)),
            top_n,
            modulo_n,
        };

        for (gene_id, gene_sequence) in genes {
            let sequence = gene_sequence.as_bytes().to_vec();
            if sequence.len() < kmer_size {
                panic!("Prefiltering: gene {} sequence length is smaller than specified kmer_size. Reduce kmer_size value.", gene_id);
            }
            for target_pos in 0..sequence.len() - kmer_size + 1 {
                let kmer: Kmer = sequence[target_pos..target_pos + kmer_size].to_vec();
                let vec = prefiltering
                    .kmer_genes_lookup
                    .entry(kmer)
                    .or_insert(Vec::new());
                vec.push(KmerGeneIndexEntry {
                    gene_id: gene_id.to_string(),
                    position: target_pos,
                });
            }
        }
        prefiltering
    }

    pub fn calculate_top_matches_with_rev_comp(
        &self,
        query: String,
    ) -> PyResult<PrefilteringResult> {
        if query.len() < self.kmer_size {
            return Err(PyValueError::new_err(
                "Query sequence shorter than kmer_size - reduce kmer_size.",
            ));
        }
        let query_bytes = query.as_bytes().to_vec();
        let top_matches = self
            .find_top_matches(&query_bytes)
            .map(|(gene_id, coverage)| GeneMatch {
                gene_id,
                rev_comp: false,
                coverage,
            });

        let rev_comp_query_bytes = dna::revcomp(&query_bytes);
        let rev_comp_top_matches =
            self.find_top_matches(&rev_comp_query_bytes)
                .map(|(gene_id, coverage)| GeneMatch {
                    gene_id,
                    rev_comp: true,
                    coverage,
                });

        let top_matches = top_matches
            .chain(rev_comp_top_matches)
            .sorted()
            .take(self.top_n)
            .collect_vec();
        let rev_comp_query = std::str::from_utf8(&rev_comp_query_bytes)
            .unwrap()
            .to_string();
        Ok(PrefilteringResult {
            query,
            rev_comp_query,
            top_matches,
        })
    }

    pub fn calculate_top_matches(&self, query: String) -> PyResult<PrefilteringResult> {
        if query.len() < self.kmer_size {
            return Err(PyValueError::new_err(
                "Query sequence shorter than kmer_size - reduce kmer_size.",
            ));
        }
        let query_byte = query.as_bytes().to_vec();
        let top_matches = self
            .find_top_matches(&query_byte)
            .map(|(gene_id, coverage)| GeneMatch {
                gene_id,
                rev_comp: false,
                coverage,
            })
            .sorted()
            .take(self.top_n)
            .collect_vec();

        let rev_comp_query = "-".to_string();
        Ok(PrefilteringResult {
            query,
            rev_comp_query,
            top_matches,
        })
    }
}

impl Prefiltering {
    fn calculate_coverage(&self, matches: Vec<&KmerMatch>) -> Coverage {
        let ksize = self.kmer_size as i32;
        let mut coverage = 0;
        if matches.is_empty() {
            return coverage;
        }
        for i in 0..matches.len() - 1 {
            let target_pos = matches[i].target_position as i32;
            let next_target_pos = matches[i + 1].target_position as i32;

            let overlap = target_pos + ksize - next_target_pos;
            match overlap < 0 {
                true => coverage += ksize,
                false => coverage += ksize - overlap,
            }
        }
        // Add last kmer
        coverage += ksize;
        coverage
    }

    fn filter_diagonal<'a>(
        &self,
        matches: Vec<&'a KmerMatch>,
    ) -> (Vec<&'a KmerMatch>, Vec<&'a KmerMatch>) {
        let mut diagonal: Vec<&KmerMatch> = Vec::new();
        let mut other_matches: Vec<&KmerMatch> = Vec::new();
        let mut prev_kmer_match = matches[0];
        diagonal.push(prev_kmer_match);

        for &kmer_match in &matches[1..] {
            let target_dist =
                kmer_match.target_position as i32 - prev_kmer_match.target_position as i32;
            let query_dist =
                kmer_match.query_position as i32 - prev_kmer_match.query_position as i32;

            match target_dist.abs_diff(query_dist) as i32 <= self.distance_threshold {
                true => {
                    diagonal.push(kmer_match);
                    prev_kmer_match = kmer_match;
                }
                false => other_matches.push(kmer_match),
            };
        }
        (diagonal, other_matches)
    }

    fn get_coverage_recursive(&self, matches: Vec<&KmerMatch>) -> Coverage {
        if matches.is_empty() {
            return 0;
        }
        let (diagonal, other_matches) = self.filter_diagonal(matches.clone());

        let current_coverage = self.calculate_coverage(diagonal);
        if current_coverage > other_matches.len() as i32 * self.kmer_size as i32 {
            return current_coverage;
        }
        let filtered_out_coverage = self.get_coverage_recursive(other_matches);
        if current_coverage < filtered_out_coverage {
            filtered_out_coverage
        } else {
            current_coverage
        }
    }

    fn find_top_matches(
        &self,
        query: &ByteSequence,
    ) -> Box<dyn Iterator<Item = (GeneId, Coverage)> + '_> {
        let mut matches: AHashMap<GeneId, Vec<KmerMatch>> =
            AHashMap::with_hasher(RandomState::with_seed(RANDOM_SEED));
        for query_pos in (0..query.len() - self.kmer_size + 1).step_by(self.modulo_n) {
            let kmer = &query[query_pos..query_pos + self.kmer_size];
            let genes: Option<&Vec<KmerGeneIndexEntry>> = self.kmer_genes_lookup.get(kmer);
            match genes {
                None => continue,
                Some(genes) => {
                    for kmer_gene_index_entry in genes {
                        let match_query_positions = matches
                            .entry(kmer_gene_index_entry.gene_id.clone())
                            .or_default();
                        match_query_positions.push(KmerMatch {
                            query_position: query_pos,
                            target_position: kmer_gene_index_entry.position,
                        });
                    }
                }
            }
        }

        Box::new(matches.into_iter().map(|(gene, mut matches)| {
            matches.sort();
            (
                gene,
                self.get_coverage_recursive(matches.iter().collect_vec()),
            )
        }))
    }
}
