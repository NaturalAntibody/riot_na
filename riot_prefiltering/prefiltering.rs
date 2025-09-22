use pyo3::{exceptions::PyValueError, prelude::*};

use crate::model::{
    ByteSequence, Coverage, GeneId, GeneMatch, GeneSegment, Kmer, KmerGeneIndexEntry, KmerMatch,
    PrefilteringResult, PrefilteringSegmentResult, SegmentMatch,
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
    min_coverage: usize,
    min_segment_length: usize,
    kmer_genes_lookup: AHashMap<Kmer, Vec<KmerGeneIndexEntry>>,
}

#[pymethods]
impl Prefiltering {
    #[new]
    #[pyo3(signature = (genes, kmer_size, distance_threshold, top_n, modulo_n, min_segment_length=None, min_coverage=None))]
    pub fn new(
        genes: HashMap<String, String>,
        kmer_size: usize,
        distance_threshold: i32,
        top_n: usize,
        modulo_n: usize,
        min_segment_length: Option<usize>,
        min_coverage: Option<usize>,
    ) -> Self {
        let resolved_min_segment_length = min_segment_length.unwrap_or(0);
        let resolved_min_coverage = min_coverage.unwrap_or(20);

        let mut prefiltering: Prefiltering = Self {
            kmer_size,
            distance_threshold,
            kmer_genes_lookup: AHashMap::with_hasher(RandomState::with_seed(RANDOM_SEED)),
            top_n,
            modulo_n,
            min_coverage: resolved_min_coverage,
            min_segment_length: resolved_min_segment_length,
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

    /// Validate query length and return bytes
    fn validate_and_convert_query(&self, query: &str) -> PyResult<Vec<u8>> {
        if query.len() < self.kmer_size {
            return Err(PyValueError::new_err(
                "Query sequence shorter than kmer_size - reduce kmer_size.",
            ));
        }
        Ok(query.as_bytes().to_vec())
    }

    pub fn calculate_top_matches_with_rev_comp(
        &self,
        query: String,
    ) -> PyResult<PrefilteringResult> {
        let query_bytes = self.validate_and_convert_query(&query)?;

        let forward_matches = self
            .find_top_matches(&query_bytes)
            .map(|(gene_id, coverage)| GeneMatch {
                gene_id,
                rev_comp: false,
                coverage,
            });

        let rev_comp_query_bytes = dna::revcomp(&query_bytes);
        let reverse_matches = self
            .find_top_matches(&rev_comp_query_bytes)
            .map(|(gene_id, coverage)| GeneMatch {
                gene_id,
                rev_comp: true,
                coverage,
            });

        let top_matches = forward_matches
            .chain(reverse_matches)
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
        let query_bytes = self.validate_and_convert_query(&query)?;

        let top_matches = self
            .find_top_matches(&query_bytes)
            .map(|(gene_id, coverage)| GeneMatch {
                gene_id,
                rev_comp: false,
                coverage,
            })
            .sorted()
            .take(self.top_n)
            .collect_vec();

        Ok(PrefilteringResult {
            query,
            rev_comp_query: "-".to_string(),
            top_matches,
        })
    }

    /// Calculate segment-centric matches (forward orientation only)
    pub fn calculate_segment_matches(&self, query: String) -> PyResult<PrefilteringSegmentResult> {
        let query_bytes = self.validate_and_convert_query(&query)?;
        let segments = self.find_non_overlapping_segments(&query_bytes, false);

        Ok(PrefilteringSegmentResult {
            query,
            rev_comp_query: "-".to_string(),
            segments,
        })
    }

    /// Calculate segment-centric matches (both orientations)
    pub fn calculate_segment_matches_with_rev_comp(&self, query: String) -> PyResult<PrefilteringSegmentResult> {
        let query_bytes = self.validate_and_convert_query(&query)?;
        let rev_comp_query_bytes = dna::revcomp(&query_bytes);

        // Get segments from both orientations
        let forward_segments = self.find_non_overlapping_segments(&query_bytes, false);
        let reverse_segments = self.find_non_overlapping_segments(&rev_comp_query_bytes, true);

        // Combine and select best non-overlapping segments
        let combined_segments = self.select_best_non_overlapping_segments(forward_segments, reverse_segments);

        let rev_comp_query = std::str::from_utf8(&rev_comp_query_bytes)
            .unwrap()
            .to_string();

        Ok(PrefilteringSegmentResult {
            query,
            rev_comp_query,
            segments: combined_segments,
        })
    }
}

impl Prefiltering {
    pub(crate) fn calculate_coverage(&self, matches: Vec<&KmerMatch>) -> Coverage {
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

    /// Finds all separate diagonal regions from a set of matches
    fn find_all_diagonals<'a>(
        &self,
        mut matches: Vec<&'a KmerMatch>,
    ) -> Vec<Vec<&'a KmerMatch>> {
        if matches.is_empty() {
            return vec![];
        }

        let mut diagonals: Vec<Vec<&KmerMatch>> = vec![];

        while !matches.is_empty() {
            let (diagonal, remaining) = self.filter_diagonal(matches);
            diagonals.push(diagonal);
            matches = remaining;
        }

        diagonals
    }

    fn filter_diagonal<'a>(
        &self,
        matches: Vec<&'a KmerMatch>,
    ) -> (Vec<&'a KmerMatch>, Vec<&'a KmerMatch>) {
        if matches.is_empty() {
            return (vec![], vec![]);
        }

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

    /// Find all kmer matches for a query sequence
    fn find_all_kmer_matches(&self, query: &ByteSequence) -> AHashMap<GeneId, Vec<KmerMatch>> {
        let mut matches: AHashMap<GeneId, Vec<KmerMatch>> =
            AHashMap::with_hasher(RandomState::with_seed(RANDOM_SEED));

        if query.len() < self.kmer_size {
            return matches;
        }

        for query_pos in (0..query.len() - self.kmer_size + 1).step_by(self.modulo_n) {
            let kmer = &query[query_pos..query_pos + self.kmer_size];
            if let Some(genes) = self.kmer_genes_lookup.get(kmer) {
                for kmer_gene_index_entry in genes {
                    matches
                        .entry(kmer_gene_index_entry.gene_id.clone())
                        .or_default()
                        .push(KmerMatch {
                            query_position: query_pos,
                            target_position: kmer_gene_index_entry.position,
                        });
                }
            }
        }

        matches
    }

    pub(crate) fn find_top_matches(
        &self,
        query: &ByteSequence,
    ) -> Box<dyn Iterator<Item = (GeneId, Coverage)> + '_> {
        let matches = self.find_all_kmer_matches(query);

        Box::new(matches.into_iter().map(|(gene, mut matches)| {
            matches.sort();
            (
                gene,
                self.get_coverage_recursive(matches.iter().collect_vec()),
            )
        }))
    }

    /// Check if a segment meets the minimum coverage requirement
    fn has_enough_coverage(&self, segment: &GeneSegment) -> bool {
        segment.coverage >= self.min_coverage as i32
    }

    /// Get all segments with both filtered and unfiltered versions for efficient processing
    fn populate_gene_segments(&self, query: &ByteSequence) -> Vec<GeneSegment> {
        let matches = self.find_all_kmer_matches(query);

        let mut filtered_gene_segments = Vec::new();

        for (gene, mut gene_matches) in matches {
            gene_matches.sort();
            let match_refs = gene_matches.iter().collect_vec();

            let gene_segments = self.convert_matches_into_gene_segments(match_refs.clone(), &gene);
            filtered_gene_segments.extend(gene_segments);
        }

        filtered_gene_segments
    }

    /// Convert gene segments to segment matches
    fn populate_single_gene_segments(&self, gene_segments: &Vec<GeneSegment>, rev_comp: bool) -> Vec<SegmentMatch> {
        gene_segments
            .iter()
            .map(|gene_segment| {
                let gene_match = GeneMatch {
                    gene_id: gene_segment.gene_id.clone(),
                    rev_comp,
                    coverage: gene_segment.coverage,
                };


                SegmentMatch {
                    query_start: gene_segment.start_query,
                    query_end: gene_segment.end_query,
                    coverage: gene_segment.coverage,
                    match_count: gene_segment.match_count,
                    matching_genes: vec![gene_match],
                    min_target_start: gene_segment.start_target,
                    max_target_start: gene_segment.start_target,
                    segment_start: gene_segment.start_query.saturating_sub(gene_segment.start_target),
                }
            })
            .collect()
    }

    /// Creates a GeneSegment from a set of diagonal matches
    fn create_gene_segment_from_diagonal(&self, diagonal: &[&KmerMatch], gene_id: &GeneId) -> GeneSegment {
        if diagonal.is_empty() {
            return GeneSegment {
                start_target: 0,
                end_target: 0,
                start_query: 0,
                end_query: 0,
                coverage: 0,
                match_count: 0,
                gene_id: "".to_string(),
            };
        }

        let start_target = diagonal.iter().map(|m| m.target_position).min().unwrap();
        let end_target = diagonal.iter().map(|m| m.target_position).max().unwrap() + self.kmer_size;
        let start_query = diagonal.iter().map(|m| m.query_position).min().unwrap();
        let end_query = diagonal.iter().map(|m| m.query_position).max().unwrap() + self.kmer_size;

        // Use the same coverage calculation method as the original API
        let coverage = self.calculate_coverage(diagonal.to_vec());
        let match_count = diagonal.len();

        GeneSegment {
            start_target,
            end_target,
            start_query,
            end_query,
            coverage,
            match_count,
            gene_id: gene_id.clone(),
        }
    }

    /// Converts diagonal matches to segments and filters by minimum coverage
    fn convert_matches_into_gene_segments(&self, matches: Vec<&KmerMatch>, gene_id: &GeneId) -> Vec<GeneSegment> {
        if matches.is_empty() {
            return vec![];
        }

        // Find all separate diagonal regions
        let diagonals = self.find_all_diagonals(matches.clone());

        // Convert each diagonal to a segment and filter by coverage in one pass
        let mut gene_segments: Vec<GeneSegment> = diagonals
            .iter()
            .filter(|diagonal| !diagonal.is_empty())
            .map(|diagonal| self.create_gene_segment_from_diagonal(diagonal, gene_id))
            .filter(|segment| self.has_enough_coverage(segment))
            .collect();

        // Prevent overlapping segments inside a gene
        // by choosing only the best one as merging them isn't clear enough
        // Sort by query start position to process segments in order
        gene_segments.sort_by(|a, b| a.start_query.cmp(&b.start_query));

        if gene_segments.is_empty() {
            return vec![];
        }

        let mut final_segments: Vec<GeneSegment> = Vec::new();
        let mut segments_iter = gene_segments.into_iter();
        let mut current_best_segment = segments_iter.next().unwrap();

        for next_segment in segments_iter {
            if next_segment.start_query < current_best_segment.end_query {
                // Overlap detected, choose the one with higher coverage
                if next_segment.coverage > current_best_segment.coverage {
                    current_best_segment = next_segment;
                }
            } else {
                // No overlap, finalize the current best segment and start a new one
                final_segments.push(current_best_segment);
                current_best_segment = next_segment;
            }
        }
        // Add the last processed segment
        final_segments.push(current_best_segment);

        final_segments
    }

    fn filter_segments_by_length(&self, segments: Vec<SegmentMatch>) -> Vec<SegmentMatch> {
        segments.into_iter().filter(|segment| self.has_enough_length_segment(segment)).collect()
    }

    fn has_enough_length_segment(&self, segment: &SegmentMatch) -> bool {
        (segment.query_end - segment.query_start) >= self.min_segment_length as usize
    }

    /// Find non-overlapping segments with gene mappings
    pub(crate) fn find_non_overlapping_segments(&self, query: &ByteSequence, rev_comp: bool) -> Vec<SegmentMatch> {
        // Get filtered gene segments
        let gene_segments = self.populate_gene_segments(query);
        let single_gene_segments = self.populate_single_gene_segments(&gene_segments, rev_comp);

        let merged_segments = self.merge_overlapping_segments(single_gene_segments);
        self.filter_segments_by_length(merged_segments)
    }

    /// Combine segments from forward and reverse complement orientations
    fn select_best_non_overlapping_segments(
        &self,
        forward_segments: Vec<SegmentMatch>,
        reverse_segments: Vec<SegmentMatch>,
    ) -> Vec<SegmentMatch> {
        // Combine all segments
        let mut all_segments = forward_segments;
        all_segments.extend(reverse_segments);

        // Use simple segment extension to capture wider coverage
        let merged_segments = self.merge_overlapping_segments(all_segments);
        self.filter_segments_by_length(merged_segments)
    }


    /// Simple segment extension - merge overlapping segments by extending boundaries
    fn merge_overlapping_segments(&self, segments: Vec<SegmentMatch>) -> Vec<SegmentMatch> {
        if segments.is_empty() {
            return vec![];
        }

        // Sort segments by segment start position (not query_start)
        let mut sorted_segments = segments;
        sorted_segments.sort_by(|a, b| a.segment_start.cmp(&b.segment_start));

        let mut extended_segments = Vec::new();
        let mut segments_iter = sorted_segments.into_iter();
        let mut current_segment = segments_iter.next().unwrap();

        for next_segment in segments_iter {
            // Check if segments overlap in the query sequence
            // Segments overlap if: next_segment.query_start <= current_segment.query_end
            if next_segment.segment_start <= current_segment.query_end {
                // Extend current segment to include the next one
                current_segment = self.extend_segment_boundaries(current_segment, next_segment);
            } else {
                // No overlap - keep current and start with next
                extended_segments.push(current_segment);
                current_segment = next_segment;
            }
        }
        // Don't forget the last segment
        extended_segments.push(current_segment);

        // Filter matching_genes to top n by coverage for each extended segment
        for segment in &mut extended_segments {
            // Keep only the top_n genes
            if segment.matching_genes.len() > self.top_n {
                segment.matching_genes.sort();
                segment.matching_genes.truncate(self.top_n);
                segment.recalculate_bounds_and_coverage();
            }
        }

        extended_segments
    }

    /// Extend segment range to include another segment's coverage
    fn extend_segment_boundaries(&self, seg1: SegmentMatch, seg2: SegmentMatch) -> SegmentMatch {
        // Extend the query boundaries to cover both segments
        let query_start = seg1.query_start.min(seg2.query_start);
        let query_end = seg1.query_end.max(seg2.query_end);
        let coverage = seg1.coverage.max(seg2.coverage);
        let min_target_start = seg1.min_target_start.min(seg2.min_target_start);
        let max_target_start = seg1.max_target_start.max(seg2.max_target_start);
        let segment_start = seg1.segment_start.min(seg2.segment_start);

        // Combine all genes from both segments
        let mut all_genes = seg1.matching_genes;
        all_genes.extend(seg2.matching_genes);
        all_genes.dedup();
        all_genes.sort();

        SegmentMatch {
            query_start: query_start,
            query_end: query_end,
            coverage: coverage,
            match_count: all_genes.len(),
            matching_genes: all_genes,
            min_target_start: min_target_start,
            max_target_start: max_target_start,
            segment_start: segment_start,
        }
    }
}


