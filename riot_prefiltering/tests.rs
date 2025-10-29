#[cfg(test)]
mod prefiltering_tests {
    use bio::alphabets::dna;
    use crate::prefiltering::Prefiltering;
    use crate::model::KmerMatch;
    use std::collections::HashMap;

    #[test]
    fn test_core_kmer_matching() {
        // Test basic kmer matching functionality
        let mut genes = HashMap::new();
        genes.insert("TEST_GENE".to_string(), "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_string());

        let prefiltering = Prefiltering::new(
            genes,
            6,  // kmer_size
            1,  // distance_threshold
            5,  // top_n
            1,  // modulo_n
            None, // min_segment_length
            None, // min_coverage
        );

        let query_bytes = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".as_bytes().to_vec();
        let matches: Vec<(String, i32)> = prefiltering.find_top_matches(&query_bytes).collect();

        assert!(!matches.is_empty());
        assert!(matches.iter().any(|(gene_id, _)| gene_id == "TEST_GENE"));
    }

    #[test]
    fn test_amino_acids_segment_matching() {
        // Test basic kmer matching functionality
        let mut genes = HashMap::new();
        genes.insert("TEST_GENE".to_string(), "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYAMHWVRQAPGQRLEWMGWINAGNGNTKYSQKFQGRVTITRDTSASTAYMELSSLRSEDTAVYYCAR".to_string());

        let prefiltering = Prefiltering::new(
            genes,
            3,  // kmer_size
            3,  // distance_threshold
            5,  // top_n
            1,  // modulo_n
            Some(60), // min_segment_length
            Some(20), // min_coverage
        );

        let query_bytes = "EIQLQQSGPELGKPGASVKVSCRASGFSFADYYIYWVKQSHGKSLELIGYIDPFNGGDTYNQIFKGKATLTVDKSSSTAFMYLNSLTSEDSAVYYCAAFRNPSFDFWGQGTTLTVSS".as_bytes().to_vec();
        let matches = prefiltering.find_non_overlapping_segments(&query_bytes, false);

        assert!(!matches.is_empty());
    }

    #[test]
    fn test_amino_acids_multiple_segment_matching() {
        // Test basic kmer matching functionality
        let mut genes = HashMap::new();

        // 'human|igk|IGKV1-17*03'
        genes.insert("IGHV3-23*01".to_string(), "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK".to_string());
        genes.insert("IGHV3-23*03".to_string(), "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSVIYSGGSSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK".to_string());

        genes.insert("IGKV1D-39*01".to_string(), "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTP".to_string());
        genes.insert("IGKV1-12*01".to_string(), "DIQMTQSPSSVSASVGDRVTITCRASQGISSWLAWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQANSFP".to_string());
        genes.insert("IGKV1-9*01".to_string(), "DIQLTQSPSFLSASVGDRVTITCRASQGISSYLAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTEFTLTISSLQPEDFATYYCQQLNSYP".to_string());
        genes.insert("IGKV1-9*02".to_string(), "DIQLTQSPSFLSASVGDRVTITCWASQGISSYLAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTEFTLTISSLQPEDFATYYCQQLNSYP".to_string());
        genes.insert("IGKV1D-16*01".to_string(), "DIQMTQSPSSLSASVGDRVTITCRASQGISSWLAWYQQKPEKAPKSLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQYNSYP".to_string());
        genes.insert("IGKV1D-16*02".to_string(), "DIQMTQSPSSLSASVGDRVTITCRARQGISSWLAWYQQKPEKAPKSLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQYNSYP".to_string());
        genes.insert("IGKV1-27*01".to_string(), "DIQMTQSPSSLSASVGDRVTITCRASQGISNYLAWYQQKPGKVPKLLIYAASTLQSGVPSRFSGSGSGTDFTLTISSLQPEDVATYYCQKYNSAP".to_string());
        genes.insert("IGKV1-6*01".to_string(), "AIQMTQSPSSLSASVGDRVTITCRASQGIRNDLGWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCLQDYNYP".to_string());
        genes.insert("IGKV1-16*02".to_string(), "DIQMTQSPSSLSASVGDRVTITCRASQGISNYLAWFQQKPGKAPKSLIYAASSLQSGVPSKFSGSGSGTDFTLTISSLQPEDFATYYCQQYNSYP".to_string());
        genes.insert("IGKV1-17*01".to_string(), "DIQMTQSPSSLSASVGDRVTITCRASQGIRNDLGWYQQKPGKAPKRLIYAASSLQSGVPSRFSGSGSGTEFTLTISSLQPEDFATYYCLQHNSYP".to_string());
        genes.insert("IGKV1-5*01".to_string(), "DIQMTQSPSTLSASVGDRVTITCRASQSISSWLAWYQQKPGKAPKLLIYDASSLESGVPSRFSGSGSGTEFTLTISSLQPDDFATYYCQQYNSYS".to_string());
        genes.insert("IGKV1-5*03".to_string(), "DIQMTQSPSTLSASVGDRVTITCRASQSISSWLAWYQQKPGKAPKLLIYKASSLESGVPSRFSGSGSGTEFTLTISSLQPDDFATYYCQQYNSYS".to_string());
        genes.insert("IGKV1-8*03".to_string(), "AIRMTQSPSSLSASTGDRVTITCRASQGISSYLAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTDFTLTISCLQSEDFATYYCQQYYSYP".to_string());
        genes.insert("IGKV1-NL1*01".to_string(), "DIQMTQSPSSLSASVGDRVTITCRASQGISNSLAWYQQKPGKAPKLLLYAASRLESGVPSRFSGSGSGTDYTLTISSLQPEDFATYYCQQYYSTP".to_string());

        genes.insert("IGKV1-8*01".to_string(), "AIRMTQSPSSFSASTGDRVTITCRASQGISSYLAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTDFTLTISCLQSEDFATYYCQQYYSYP".to_string());
        genes.insert("IGKV1-8*02".to_string(), "AIRITQSPSSLSASTGDRVTITCRASQGISSYLAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTDFTLTISCLQSEDFATYYCQQYYSYP".to_string());
        genes.insert("IGKV1D-13*02".to_string(), "AIQLTQSPSSLSASVGDRVTITCRASQGISSALAWYQQKPGKAPKLLIYDASSLESGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQFNSYP".to_string());
        genes.insert("IGKV1D-33*01".to_string(), "DIQMTQSPSSLSASVGDRVTITCQASQDISNYLNWYQQKPGKAPKLLIYDASNLETGVPSRFSGSGSGTDFTFTISSLQPEDIATYYCQQYDNLP".to_string());
        genes.insert("IGKV1D-43*01".to_string(), "AIRMTQSPFSLSASVGDRVTITCWASQGISSYLAWYQQKPAKAPKLFIYYASSLQSGVPSRFSGSGSGTDYTLTISSLQPEDFATYYCQQYYSTP".to_string());

        genes.insert("IGKV1-17*03".to_string(), "DIQMTQSPSAMSASVGDRVTITCRASQGISNYLAWFQQKPGKVPKRLIYAASSLQSGVPSRFSGSGSGTEFTLTISSLQPEDFATYYCLQHNSYP".to_string());

        let prefiltering = Prefiltering::new(
            genes,
            3,  // kmer_size
            4,  // distance_threshold
            20,  // top_n
            1,  // modulo_n
            Some(60), // min_segment_length
            Some(20), // min_coverage
        );

        let query_bytes = "MPSSAVGVLGEAWYSLGGPDSSCAASGFTFSSYAMSWVRQAPGKGLEWVSSIANKGHETRYVDSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKYAGTFDYWGQGTLVTVSSGGGGSGGGGSGGGGSTDIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASMLQSGVPSRFSGSGSGTDFTLTISSLQPEYFATYYCQQARSWPPTFGQGDQGGNQTGRPHIIIAITGATHHHHHHGAAEQKLISEEDLNGA".as_bytes().to_vec();
        let matches = prefiltering.find_non_overlapping_segments(&query_bytes, false);

        println!("matches ({}):", matches.len());
        for (i, seg) in matches.iter().enumerate() {
            println!("  [{}] {:?}", i, seg);
        }

        assert!(!matches.is_empty());
    }

    #[test]
    fn test_short_kmer_matching() {
        // Test basic kmer matching functionality
        let mut genes = HashMap::new();
        genes.insert("TEST_GENE".to_string(), "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_string());

        let prefiltering = Prefiltering::new(
            genes,
            6,  // kmer_size
            1,  // distance_threshold
            5,  // top_n
            1,  // modulo_n
            None, // min_segment_length
            None, // min_coverage
        );

        let query_bytes = "ACGT".as_bytes().to_vec();
        let matches: Vec<(String, i32)> = prefiltering.find_top_matches(&query_bytes).collect();

        assert!(matches.is_empty());
    }

    #[test]
    fn test_coverage_calculation() {
        // Test coverage calculation with known inputs
        let mut genes = HashMap::new();
        genes.insert("GENE1".to_string(), "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_string());

        let prefiltering = Prefiltering::new(
            genes,
            6,  // kmer_size
            1,  // distance_threshold
            5,  // top_n
            1,  // modulo_n
            None,
            None,
        );

        // Test coverage calculation with known kmer matches
        let matches = vec![
            &KmerMatch { target_position: 0, query_position: 0 },
            &KmerMatch { target_position: 1, query_position: 1 },
            &KmerMatch { target_position: 2, query_position: 2 },
        ];

        let coverage = prefiltering.calculate_coverage(matches);
        assert!(coverage > 0);
    }

    #[test]
    fn test_empty_and_edge_inputs() {
        // Test edge cases with empty or minimal inputs
        let mut genes = HashMap::new();
        genes.insert("MINIMAL_GENE".to_string(), "ACGTACGT".to_string()); // Minimal gene sequence

        let prefiltering = Prefiltering::new(
            genes,
            6,  // kmer_size
            1,  // distance_threshold
            5,  // top_n
            1,  // modulo_n
            None,
            None,
        );

        // Test with minimal query
        let query_bytes = "ACGTAC".as_bytes().to_vec();
        let matches: Vec<(String, i32)> = prefiltering.find_top_matches(&query_bytes).collect();

        // Should handle edge case gracefully (may be empty or have results)
        assert!(matches.len() <= 5); // Should respect top_n limit
    }

    #[test]
    fn test_reverse_complement_functionality() {
        // Test reverse complement handling in core functionality
        let mut genes = HashMap::new();
        let gene_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        genes.insert("FORWARD_GENE".to_string(), gene_seq.to_string());
        genes.insert("REVERSE_GENE".to_string(), String::from_utf8(dna::revcomp(gene_seq.as_bytes().to_vec())).unwrap());

        let prefiltering = Prefiltering::new(
            genes,
            6,  // kmer_size
            1,  // distance_threshold
            10,  // top_n
            1,  // modulo_n
            None,
            None,
        );

        let query_bytes = gene_seq.as_bytes().to_vec();
        let forward_matches: Vec<(String, i32)> = prefiltering.find_top_matches(&query_bytes).collect();

        let rev_comp_query_bytes = dna::revcomp(&query_bytes);
        let reverse_matches: Vec<(String, i32)> = prefiltering.find_top_matches(&rev_comp_query_bytes).collect();

        assert!(!forward_matches.is_empty());
        assert!(!reverse_matches.is_empty());
    }

    #[test]
    fn test_segment_overlap_detection() {
        // Test core segment overlap detection logic
        let mut genes = HashMap::new();
        genes.insert("GENE1".to_string(), "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_string());
        genes.insert("GENE2".to_string(), "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA".to_string());

        let prefiltering = Prefiltering::new(
            genes,
            6,  // kmer_size
            1,  // distance_threshold
            5,  // top_n
            1,  // modulo_n
            Some(10), // min_segment_length
            None,
        );

        let query = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA";
        let query_bytes = query.as_bytes().to_vec();
        let segments = prefiltering.find_non_overlapping_segments(&query_bytes, false);

        // Should find non-overlapping segments
        assert!(!segments.is_empty());

        // Check that segments don't overlap (basic validation)
        let mut segments_sorted = segments.clone();
        segments_sorted.sort_by_key(|s| s.query_start);
        for i in 0..segments_sorted.len().saturating_sub(1) {
            assert!(segments_sorted[i].query_end <= segments_sorted[i + 1].query_start,
                   "Segments should not overlap");
        }
    }

    #[test]
    fn test_segment_overlap_detection_with_min_segment_length() {
        // Test core segment overlap detection logic with min_segment_length
        let mut genes = HashMap::new();
        let gene_seq = "GAAGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGCAGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTGATGATTATGCCATGCACTGGGTCCGGCAAGCTCCAGGGAAGGGCCTGGAGTGGGTCTCAGGTATTAGTTGGAATAGTGGTAGCATAGGCTATGCGGACTCTGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCCCTGTATCTGCAAATGAACAGTCTGAGAGCTGAGGACACGGCCTTGTATTACTGTGCAAAAGATA";
        genes.insert("GENE1".to_string(), gene_seq.to_string());
        //genes.insert("GENE2".to_string(), String::from_utf8(dna::revcomp(gene_seq.as_bytes().to_vec())).unwrap());

        let prefiltering = Prefiltering::new(
            genes,
            9,  // kmer_size
            13,  // distance_threshold
            5,  // top_n
            2,  // modulo_n
            Some(60), // min_segment_length
            None,
        );

        let query = "CTTCCGATCTATTTCCCTTAGACTCTCCTGTGCAGCGTCTGGATTCAGAGACAACGCCAAGAACTCCCTGTATCTGGAAATGAACAGTCTGAGACCTGAGGACACGGCCCTCTATTACTGTGTGAAAGATCGGGGTGGTGCGGGGAATTCGCCCACGTACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAGCCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCGTACGTCCGAGATCGG";
        let query_bytes = query.as_bytes().to_vec();
        let segments = prefiltering.find_non_overlapping_segments(&query_bytes, true);

        // Should find non-overlapping segments
        assert!(!segments.is_empty());

        // Check that segments don't overlap (basic validation)
        let mut segments_sorted = segments.clone();
        assert_eq!(segments_sorted.len(), 1);
        segments_sorted.sort_by_key(|s| s.query_start);
        for i in 0..segments_sorted.len().saturating_sub(1) {
            assert!(segments_sorted[i].query_end <= segments_sorted[i + 1].query_start,
                   "Segments should not overlap");
        }
    }

    #[test]
    fn test_kmer_index_construction() {
        // Test that the kmer index is built correctly
        let mut genes = HashMap::new();
        genes.insert("SHORT_GENE".to_string(), "ACGTACGT".to_string());

        // This should work without panicking
        let prefiltering = Prefiltering::new(
            genes,
            6,  // kmer_size smaller than gene
            1,  // distance_threshold
            5,  // top_n
            1,  // modulo_n
            None,
            None,
        );

        // Basic functionality test
        let query_bytes = "ACGTAC".as_bytes().to_vec();
        let matches: Vec<(String, i32)> = prefiltering.find_top_matches(&query_bytes).collect();

        // Should handle this correctly
        assert!(matches.len() <= 5);
    }

    #[test]
    fn test_modulo_stepping() {
        // Test that modulo stepping works correctly
        let mut genes = HashMap::new();
        genes.insert("TEST_GENE".to_string(), "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_string());

        let prefiltering_mod1 = Prefiltering::new(
            genes.clone(),
            6,  // kmer_size
            1,  // distance_threshold
            10, // top_n
            1,  // modulo_n = 1 (check every position)
            None,
            None,
        );

        let prefiltering_mod3 = Prefiltering::new(
            genes,
            6,  // kmer_size
            1,  // distance_threshold
            10, // top_n
            3,  // modulo_n = 3 (check every 3rd position)
            None,
            None,
        );

        let query_bytes = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".as_bytes().to_vec();
        let matches_mod1: Vec<(String, i32)> = prefiltering_mod1.find_top_matches(&query_bytes).collect();
        let matches_mod3: Vec<(String, i32)> = prefiltering_mod3.find_top_matches(&query_bytes).collect();

        // Both should find matches, but mod1 might be more sensitive
        assert!(!matches_mod1.is_empty());
        assert!(!matches_mod3.is_empty());
    }
}