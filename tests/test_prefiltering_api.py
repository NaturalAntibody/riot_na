from typing import Dict

import pytest

import riot_na


class TestPrefilteringNucleotideAPI:
    """Test suite for the Prefiltering API with nucleotide sequences."""

    @pytest.fixture
    def basic_nucleotide_genes(self) -> Dict[str, str]:
        """Basic nucleotide gene database for testing."""
        return {
            "IGHV1-69*01": "CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTCACCAGCTATGCTATGCACTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAACCCTAACAGTGGTGGCACAAACTATGCACAGAAGTTTCAGGGCAGGGTCACCATGACCAGGGACACGTCCATCAGCACAGCCTACATGGAGCTGAGCAGGCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGA",
            "IGHV3-23*01": "GAAGTGCAGCTGCTGGAGTCTGGCGGAGGATTGGTACAGCCTGGCGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGCAGTTACTGGATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCAGCTATTAATGGTGGTAGCGGCACTTATTACGCTGACTCTGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTTTATCTGCAAATGAACAGCCTGAGAGCTGAGGACACGGCTGTGTATTATTGTGCTAAGGGTGAT",
            "IGHV4-34*01": "CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGGAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTTACCAGTTACGCTATGCACTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGCTGGATCAACGCTGGCAACGGTAACACAAAATATTCACAGAAGTTCCAGGGCAGAGTCACCATTACCAGGGACACGTCCATCAGCACAGCCTACATGGAGCTGAGCAGGCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGAT",
            "IGKV1-27*01": "GACATCCAGATGACCCAGTCTCCATCCTCCCTGTCTGCATCTGTAGGAGACAGAGTCACCATCACTTGCCGGGCAAGTCAGATCAGTAGCAGCTATTTAAATTGGTATCAGCAGAAACCAGGGAAAGCCCCTAAGCTCCTGATCTATGCTGCATCCAGTTTGCAAAGTGGGGTCCCATCAAGGTTCAGTGGCAGTGGATCTGGGACAGATTTCACTCTCACCATCAGCAGTCTGCAACCTGAAGATTTTGCAACTTATTACTGTCAGCAATCTTATAGTACTCCTCTCACTTTCGGCGGAGGGACCAAGGTGGAGATCAAACGG",
            "IGLV2-8*01": "CAGTCTGTGTTGACGCAGCCGCCCTCAGTGTCTGGGGCCCCAGGGCAGAGGGTCACCATCTCCTGCACTGGGAGCAGCTCCAACATCGGGGCAGGTTATGATGTACACTGGTACCAGCAGCTTCCAGGAACAGCCCCCAAACTCCTCATCTATGGTAACAGCAATCGGCCCTCAGGGGTCCCTGACCGATTCTCTGGCTCCAAGTCTGGCACCTCAGCCTCCCTGGCCATCAGTGGGCTCCAGTCTGAGGATGAGGCTGATTATTACTGTGCAACATGGGATGACAGCCTGAATGGTCCG",
        }

    @pytest.fixture
    def nucleotide_immunoglobulin_genes(self) -> Dict[str, str]:
        """Comprehensive nucleotide immunoglobulin gene database."""
        return {
            # Heavy chain V genes (nucleotide)
            "IGHV1-2*01": "CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTCACCAGCTATGCTATGCACTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAACCCTAACAGTGGTGGCACAAACTATGCACAGAAGTTTCAGGGCAGGGTCACCATGACCAGGGACACGTCCATCAGCACAGCCTACATGGAGCTGAGCAGGCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGA",
            "IGHV3-33*01": "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTCCAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTACGCTATGCACTGGGTCCGCCAGGCTCCAGGCAAGGGGCTGGAGTGGGTGGCAGTTATATCATATGATGGAAGTAATAAATACTATGCAGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGAGGGGGGTATA",
            "IGHV4-39*01": "CAGGTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCTTCGGAGACCCTGTCCCTCACCTGCGCTGTCTATGGTGGGTCCTTCAGTGGTTACTACTGGAGCTGGATCCGCCAGCCCCCAGGGAAGGGGCTGGAGTGGATTGGGGAAATCAATCATAGTGGAAGCACCAACTACAACCCGTCCCTCAAGAGTCGAGTCACCATATCAGTAGACACGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCGGACACGGCCGTGTATTACTGTGCGAGAGATCGGGGGAGGAA",
            # Kappa light chain V genes (nucleotide)
            "IGKV1-39*01": "GACATCCAGATGACCCAGTCTCCATCCTCCCTGTCTGCATCTGTAGGAGACAGAGTCACCATCACTTGCCGGGCAAGTCAGGGCATTAGTAGCTGGCTGGCCTGGTATCAGCAGAAACCAGGGAAAGCCCCTAAGCTCCTGATCTATGCTGCATCCAGTTTGCAAAGTGGGGTCCCATCAAGGTTCAGTGGCAGTGGATCTGGGACAGATTTCACTCTCACCATCAGCAGTCTGCAACCTGAAGATTTTGCAACTTACTACTGTCAGCAAGCTAACAGCTTCCCTTATACGTTCGGCCAAGGGACCAAGGTG",
            "IGKV2-28*01": "GATGTTGTGATGACTCAGTCTCCACTCTCCCTGCCCGTCACCCTTGGACAGCCGGCCTCCATCTCCTGCAGGTCTAGTCAAAGCCTCGTAGACAGTGATGGAAACACCTATTTGAATTGGTTGCTACAGAGGCCAGGCCAGTCTCCAAAGCGCCTAATCTATCTGGTGTCTAAACTGGACTCTGGAGTCCCTGACAGGTTCACTGGCAGTGGATCAGGGACAGATTTCACACTGAAAATCAGCAGAGTGGAGGCTGAGGATGTGGGTGTTTATTACTGCATGCAAGCTCTACAAACTCCGTACACTTTT",
            "IGKV3-15*01": "GAAATTGTGTTGACACAGTCTCCAGCCACCCTGTCTTTGTCTCCAGGGGAAAGAGCCACCCTCTCCTGCAGGGCCAGTCAGAGTGTTAGCAGCAGCTACTTAGCCTGGTACCAGCAGAAACCTGGCCAGGCTCCCAGGCTCCTCATCTATGGTGCATCCAGCAGGGCCACTGGCATCCCAGACAGGTTCAGTGGCAGTGGGTCTGGGACAGACTTCACTCTCACCATCAGCAGCCTAGAGCCTGAAGATTTTGCAGTTTATTACTGTCAGCAGTATGGTAGCTCACCTCCGTGGACGTTCGGCCAAGGG",
            # Lambda light chain V genes (nucleotide)
            "IGLV1-44*01": "CAGTCTGCCCTGACTCAGCCTGCCTCCGTGTCTGGGTCTCCTGGACAGTCGATCACCATCTCCTGCACTGGAACCAGCAGTGACGTTGGTGGTTATAACTATGTCTCCTGGTACCAACAACACCCAGGCAAAGCCCCCAAACTCATGATTTATGACGTCAGTAAGCGGCCCTCAGGGGTTTCTAATCGCTTCTCTGGCTCCAAGTCTGGCAACACGGCCTCCCTGACCATCTCTGGGCTCCAGGCTGAGGATGAGGCTGATTATTACTGCAGCTCATATACAAGCAGCAGCACTCTGGTGTTCGGCGGAGGGACCAAGCTGACCGTCCTA",
            "IGLV2-14*01": "CAGTCTGTGCTGACTCAGCCACCCTCAGCGTCTGGGACCCCCGGGCAGAGGGTCACCATCTCTTGTTCTGGAAGCAGCTCCAACATCGGAAGTAATACTGTAAACTGGTACCAGCAGCTCCCAGGAACGGCCCCCAAACTCCTCATCTATAGTAATAATCAGCGGCCCTCAGGGGTCCCTGACCGATTCTCTGGCTCCAAGTCTGGCACGTCAGCCACCCTGGGCATCACCGGACTCCAGACTGGGGACGAGGCCGATTATTACTGCGGAACATGGGATAGCAGCCTGAGTGCTGTGGTATTCGGCGGAGGGACCAAGCTGACCGTCCTA",
            # J genes (nucleotide)
            "IGHJ4*01": "ACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG",
            "IGHJ6*01": "ATTACTACTACTACTACGGTATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCAG",
            "IGKJ1*01": "TGGACGTTCGGCCAAGGGACCAAGGTGGAAATCAAACG",
            "IGLJ3*01": "TGGACGTTCGGCCAAGGGACCAAGGTGGAAATCAAACG",
        }

    @pytest.fixture
    def multi_domain_nucleotide_sequence(self) -> str:
        """A nucleotide sequence containing multiple domains for testing segment detection."""
        return (
            # Heavy chain V domain (nucleotide)
            "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTCCAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTACGCTATGCACTGGGTCCGCCAGGCTCCAGGCAAGGGGCTGGAGTGGGTGGCAGTTATATCATATGATGGAAGTAATAAATACTATGCAGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGAGGGGGGTATA"
            # Linker (nucleotide)
            "GGAGGCGGTTCAGGCGGAGGTGGCTCTGGCGGTGGCGGATCG"
            # Kappa light chain V domain (nucleotide)
            "GAAATTGTGTTGACACAGTCTCCAGCCACCCTGTCTTTGTCTCCAGGGGAAAGAGCCACCCTCTCCTGCAGGGCCAGTCAGAGTGTTAGCAGCAGCTACTTAGCCTGGTACCAGCAGAAACCTGGCCAGGCTCCCAGGCTCCTCATCTATGGTGCATCCAGCAGGGCCACTGGCATCCCAGACAGGTTCAGTGGCAGTGGGTCTGGGACAGACTTCACTCTCACCATCAGCAGCCTAGAGCCTGAAGATTTTGCAGTTTATTACTGTCAGCAGTATGGTAGCTCACCTCCGTGGACGTTCGGCCAAGGG"
            # Another linker (nucleotide)
            "GGAGGCGGTTCAGGCGGAGGTGGCTCTGGCGGTGGCGGATCG"
            # Lambda light chain V domain (nucleotide)
            "CAGTCTGCCCTGACTCAGCCTGCCTCCGTGTCTGGGTCTCCTGGACAGTCGATCACCATCTCCTGCACTGGAACCAGCAGTGACGTTGGTGGTTATAACTATGTCTCCTGGTACCAACAACACCCAGGCAAAGCCCCCAAACTCATGATTTATGACGTCAGTAAGCGGCCCTCAGGGGTTTCTAATCGCTTCTCTGGCTCCAAGTCTGGCAACACGGCCTCCCTGACCATCTCTGGGCTCCAGGCTGAGGATGAGGCTGATTATTACTGCAGCTCATATACAAGCAGCAGCACTCTGGTGTTCGGCGGAGGGACCAAGCTGACCGTCCTA"
        )

    def test_prefiltering_instantiation_basic(self, basic_nucleotide_genes):
        """Test basic Prefiltering instantiation with nucleotide genes."""
        prefiltering = riot_na.Prefiltering(
            genes=basic_nucleotide_genes, kmer_size=6, distance_threshold=1, top_n=5, modulo_n=1
        )
        assert prefiltering is not None

    def test_prefiltering_instantiation_with_min_segment_length(self, basic_nucleotide_genes):
        """Test Prefiltering instantiation with min_segment_length parameter."""
        prefiltering = riot_na.Prefiltering(
            genes=basic_nucleotide_genes, kmer_size=6, distance_threshold=1, top_n=5, modulo_n=1, min_segment_length=15
        )
        assert prefiltering is not None

    def test_calculate_top_matches_basic(self, basic_nucleotide_genes):
        """Test basic top matches calculation with nucleotide sequences."""
        prefiltering = riot_na.Prefiltering(
            genes=basic_nucleotide_genes, kmer_size=6, distance_threshold=1, top_n=3, modulo_n=1
        )

        # Use a nucleotide query that should match one of our genes
        query = "CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTCACCAGCTATGCTATGCACTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAACCCTAACAGTGGTGGCACAAACTATGCACAGAAGTTTCAGGGCAGGGTCACCATGACCAGGGACACGTCCATCAGCACAGCCTACATGGAGCTGAGCAGGCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGA"

        result = prefiltering.calculate_top_matches(query)

        # Verify result structure
        assert result is not None
        assert result.query == query
        assert result.rev_comp_query == "-"  # No reverse complement in this method
        assert isinstance(result.top_matches, list)
        assert len(result.top_matches) <= 3  # Should respect top_n limit

        # Check that at least one match was found
        if len(result.top_matches) > 0:
            match = result.top_matches[0]
            assert hasattr(match, "gene_id")
            assert hasattr(match, "coverage")
            assert hasattr(match, "rev_comp")
            assert match.rev_comp is False  # Forward orientation only

    def test_calculate_top_matches_with_rev_comp(self, basic_nucleotide_genes):
        """Test top matches calculation with reverse complement."""
        prefiltering = riot_na.Prefiltering(
            genes=basic_nucleotide_genes, kmer_size=6, distance_threshold=1, top_n=5, modulo_n=1
        )

        query = "CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTCACCAGCTATGCTATGCACTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAACCCTAACAGTGGTGGCACAAACTATGCACAGAAGTTTCAGGGCAGGGTCACCATGACCAGGGACACGTCCATCAGCACAGCCTACATGGAGCTGAGCAGGCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGA"

        result = prefiltering.calculate_top_matches_with_rev_comp(query)

        # Verify result structure
        assert result is not None
        assert result.query == query
        assert result.rev_comp_query != "-"  # Should have reverse complement
        assert len(result.rev_comp_query) == len(query)  # Same length as original
        assert isinstance(result.top_matches, list)
        assert len(result.top_matches) <= 5  # Should respect top_n limit

    def test_calculate_segment_matches_basic(self, nucleotide_immunoglobulin_genes):
        """Test basic segment matches calculation with nucleotides."""
        prefiltering = riot_na.Prefiltering(
            genes=nucleotide_immunoglobulin_genes,
            kmer_size=6,
            distance_threshold=1,
            top_n=10,
            modulo_n=1,
            min_segment_length=20,
        )

        # Use a nucleotide sequence that should match multiple genes
        query = "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTCCAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTACGCTATGCACTGGGTCCGCCAGGCTCCAGGCAAGGGGCTGGAGTGGGTGGCAGTTATATCATATGATGGAAGTAATAAATACTATGCAGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGAGGGGGGTATA"

        result = prefiltering.calculate_segment_matches(query)

        # Verify result structure
        assert result is not None
        assert result.query == query
        assert result.rev_comp_query == "-"  # No reverse complement in this method
        assert isinstance(result.segments, list)

        # Check segment structure if any found
        if len(result.segments) > 0:
            segment = result.segments[0]
            assert hasattr(segment, "query_start")
            assert hasattr(segment, "query_end")
            assert hasattr(segment, "coverage")
            assert hasattr(segment, "match_count")
            assert hasattr(segment, "matching_genes")
            assert hasattr(segment, "min_target_start")

            # Check segment bounds
            assert segment.query_start <= segment.query_end
            assert segment.query_end <= len(query)
            assert isinstance(segment.matching_genes, list)

    def test_multi_domain_segment_detection(self, nucleotide_immunoglobulin_genes, multi_domain_nucleotide_sequence):
        """Test segment detection on multi-domain nucleotide sequences."""
        prefiltering = riot_na.Prefiltering(
            genes=nucleotide_immunoglobulin_genes,
            kmer_size=6,
            distance_threshold=1,
            top_n=10,
            modulo_n=1,
            min_segment_length=15,
        )

        result = prefiltering.calculate_segment_matches_with_rev_comp(multi_domain_nucleotide_sequence)

        # Should detect multiple segments
        assert len(result.segments) > 0

        # Check that segments don't overlap (should be non-overlapping)
        segments = sorted(result.segments, key=lambda s: s.query_start)
        for i in range(len(segments) - 1):
            assert segments[i].query_end <= segments[i + 1].query_start, "Segments should not overlap"

    def test_query_length_handling(self, basic_nucleotide_genes):
        """Test behavior with nucleotide queries of different lengths."""
        prefiltering = riot_na.Prefiltering(
            genes=basic_nucleotide_genes, kmer_size=6, distance_threshold=1, top_n=5, modulo_n=1
        )

        # Test minimum valid query length (exactly kmer_size)
        min_query = "ACGTAC"  # 6 characters, exactly kmer_size
        result = prefiltering.calculate_top_matches(min_query)
        assert result is not None
        assert result.query == min_query

        # Test longer query
        long_query = "ACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        result = prefiltering.calculate_top_matches(long_query)
        assert result is not None
        assert result.query == long_query


class TestPrefilteringAminoAcidAPI:
    """Test suite for the Prefiltering API with amino acid sequences."""

    @pytest.fixture
    def basic_amino_acid_genes(self) -> Dict[str, str]:
        """Basic amino acid gene database for testing."""
        return {
            "IGHV1-69*01": "QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGRIIPILGIANYAQKFQGRVTITADKSTSTAYMELSSLRSEDTAVYYCAR",
            "IGHV3-23*01": "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK",
            "IGHV4-34*01": "QVQLVQSGGGVKKPGASVKVSCKASGYTFTSYAMHWVRQAPGQGLEWMGWINAGNGNTKYSQKFQGRVTITRDTSASTAYMELSSLRSEDTAVYYCAR",
            "IGKV1-27*01": "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPLT",
            "IGLV2-8*01": "QSALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYDVSNRPSGVSNRFSGSKSGNTASLTISGLQAEDEADYYCSSYTSSSTLV",
        }

    @pytest.fixture
    def amino_acid_immunoglobulin_genes(self) -> Dict[str, str]:
        """Comprehensive amino acid immunoglobulin gene database."""
        return {
            # Heavy chain V genes (amino acid)
            "IGHV1-2*01": "QVQLVQSGAEVKKPGSSVKVSCKASGYTFTSYAISWVRQAPGQGLEWMGRIIPILGIANYAQKFQGRVTITADKSTSTAYMELSSLRSEDTAVYYCAR",
            "IGHV1-69-2*01": "EVQLVQSGAEVKKPGATVKISCKVSGYTFTDYYMHWVQQAPGKGLEWMGLVDPEDGETIYAEKFQGRVTITADTSTDTAYMELSSLRSEDTAVYYCAT",
            "IGHV3-33*01": "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMHWVRQAPGKGLEWVAVIYHGSNKYYANPVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCAR",
            "IGHV1-3*01": "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYAMHWVRQAPGQRLEWMGWINAGNGNTKYSQKFQGRVTITRDTSASTAYMELSSLRSEDTAVYYCAR",
            "IGHV4-39*01": "QVQLQESGPGLVKPSETLSLTCAVYGGSFSGYYWSWIRQPPGKGLEWIGEINHSGSTNYPDSLKSRVISIVDTSKNQFSLKLSSVTAADTAVYYCAR",
            # Kappa light chain V genes (amino acid)
            "IGKV1-39*01": "DIQMTQSPSSLSASVGDRVTITCRASQGISSWLAWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQANSFPYT",
            "IGKV2-28*01": "DVVMTQSPLSLPVTLGQPASISCRSSQSLVHSGATYLNWFQQRPGQSPRRLLIYKVSNRDSGVPDRFTGSGSGTDFTLKISRVEAEDVGVYYCAR",
            "IGKV3-15*01": "EIVLTQSPATLSLSPGERATLSCRASQSVSSSYLAWYQQKPGQAPRLLIYGASSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQQYGSSPPWT",
            # Lambda light chain V genes (amino acid)
            "IGLV1-44*01": "QSALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYDVSNRPSGVSNRFSGSKSGNTASLTISGLQAEDEADYYCSSYTSSSTLV",
            "IGLV2-14*01": "QSALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYDVSNRPSGVSNRFSGSKSGNTASLTISGLQAEDEADYYCGTWDSSLSAVVFGGGTKLT",
            # J genes (amino acid)
            "IGHJ4*01": "YYYYGMDVWGQGTLVTVSS",
            "IGHJ6*01": "YYYYYGMDVWGQGTLVTVSS",
            "IGKJ1*01": "WTFGQGTKVEIK",
            "IGLJ3*01": "WTFGQGTKVEIK",
        }

    @pytest.fixture
    def multi_domain_amino_acid_sequence(self) -> str:
        """An amino acid sequence containing multiple domains for testing segment detection."""
        return (
            # Heavy chain V domain (amino acid)
            "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMHWVRQAPGKGLEWVAVIYHGSNKYYANPVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCAR"
            # Linker (amino acid)
            "GGGSGGGSGGGSGGS"
            # Kappa light chain V domain (amino acid)
            "EIVLTQSPATLSLSPGERATLSCRASQSVSSSYLAWYQQKPGQAPRLLIYGASSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQQYGSSPPWT"
            # Another linker (amino acid)
            "GGGSGGGSGGGSGGS"
            # Lambda light chain V domain (amino acid)
            "QSALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYDVSNRPSGVSNRFSGSKSGNTASLTISGLQAEDEADYYCSSYTSSSTLV"
        )

    def test_calculate_top_matches_basic(self, basic_amino_acid_genes):
        """Test basic top matches calculation with amino acid sequences."""
        prefiltering = riot_na.Prefiltering(
            genes=basic_amino_acid_genes,
            kmer_size=3,  # Appropriate for amino acids
            distance_threshold=1,
            top_n=3,
            modulo_n=1,
        )

        # Use an amino acid query that should match one of our genes
        query = "QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGRIIPILGIANYAQKFQGRVTITADKSTSTAYMELSSLRSEDTAVYYCAR"

        result = prefiltering.calculate_top_matches(query)

        # Verify result structure
        assert result is not None
        assert result.query == query
        assert result.rev_comp_query == "-"  # No reverse complement in this method
        assert isinstance(result.top_matches, list)
        assert len(result.top_matches) <= 3  # Should respect top_n limit

        # Check that at least one match was found
        if len(result.top_matches) > 0:
            match = result.top_matches[0]
            assert hasattr(match, "gene_id")
            assert hasattr(match, "coverage")
            assert hasattr(match, "rev_comp")
            assert match.rev_comp is False  # Forward orientation only

    def test_calculate_segment_matches_basic(self, amino_acid_immunoglobulin_genes):
        """Test basic segment matches calculation with amino acids."""
        prefiltering = riot_na.Prefiltering(
            genes=amino_acid_immunoglobulin_genes,
            kmer_size=3,  # Appropriate for amino acids
            distance_threshold=3,
            top_n=10,
            modulo_n=1,
            min_segment_length=60,  # Shorter for amino acids
        )

        # Use an amino acid sequence that should match multiple genes
        query = "EIQLQQSGPELGKPGASVKVSCRASGFSFADYYIYWVKQSHGKSLELIGYIDPFNGGDTYNQIFKGKATLTVDKSSSTAFMYLNSLTSEDSAVYYCAAFRNPSFDFWGQGTTLTVSSGGGGGGGGSGGGGGGGGGQIVLIQSPPIMSASPGEKVTLTCSASSSVSSRYLYWYQQKPGSSPKLWIYGTSNLASGVPARFSGSGSGTSFSLTISSMEAEDAASYFCHQWSSFPFTFGSGTKLEIK"

        result = prefiltering.calculate_segment_matches(query)

        # Verify result structure
        assert result is not None
        assert result.query == query
        assert result.rev_comp_query == "-"  # No reverse complement in this method
        assert isinstance(result.segments, list)
        assert len(result.segments) >= 1 and len(result.segments) <= 2

        # Check segment structure if any found
        if len(result.segments) > 0:
            segment = result.segments[0]
            assert hasattr(segment, "query_start")
            assert hasattr(segment, "query_end")
            assert hasattr(segment, "coverage")
            assert hasattr(segment, "match_count")
            assert hasattr(segment, "matching_genes")
            assert hasattr(segment, "min_target_start")

            # Check segment bounds
            assert segment.query_start <= segment.query_end
            assert segment.query_end <= len(query)
            assert isinstance(segment.matching_genes, list)

    def test_multi_domain_segment_detection(self, amino_acid_immunoglobulin_genes, multi_domain_amino_acid_sequence):
        """Test segment detection on multi-domain amino acid sequences."""
        prefiltering = riot_na.Prefiltering(
            genes=amino_acid_immunoglobulin_genes,
            kmer_size=3,  # Appropriate for amino acids
            distance_threshold=1,
            top_n=10,
            modulo_n=1,
            min_segment_length=8,
        )

        result = prefiltering.calculate_segment_matches_with_rev_comp(multi_domain_amino_acid_sequence)

        # Should detect multiple segments
        assert len(result.segments) > 0

        # Check that segments don't overlap (should be non-overlapping)
        segments = sorted(result.segments, key=lambda s: s.query_start)
        for i in range(len(segments) - 1):
            assert segments[i].query_end <= segments[i + 1].query_start, "Segments should not overlap"

    def test_therapeutic_scenarios(self, amino_acid_immunoglobulin_genes):
        """Test therapeutic scenarios with amino acid sequences."""
        prefiltering = riot_na.Prefiltering(
            genes=amino_acid_immunoglobulin_genes,
            kmer_size=3,  # Appropriate for amino acids
            distance_threshold=4,
            top_n=10,
            modulo_n=1,
            min_segment_length=60,
            min_coverage=20,
        )

        query = "EIQLQQSGPELGKPGASVKVSCRASGFSFADYYIYWVKQSHGKSLELIGYIDPFNGGDTYNQIFKGKATLTVDKSSSTAFMYLNSLTSEDSAVYYCAAFRNPSFDFWGQGTTLTVSS"
        result = prefiltering.calculate_segment_matches(query)
        assert len(result.segments) > 0
        assert result.segments[0].query_start == 5
        assert result.segments[0].query_end == 97
        assert result.segments[0].coverage == 30
        assert result.segments[0].matching_genes[0].gene_id == "IGHV1-2*01"
        assert result.segments[0].matching_genes[1].gene_id == "IGHV1-3*01"
        assert result.segments[0].matching_genes[2].gene_id == "IGHV1-69-2*01"

    def test_query_length_handling(self, basic_amino_acid_genes):
        """Test behavior with amino acid queries of different lengths."""
        prefiltering = riot_na.Prefiltering(
            genes=basic_amino_acid_genes,
            kmer_size=3,  # Appropriate for amino acids
            distance_threshold=1,
            top_n=5,
            modulo_n=1,
        )

        # Test minimum valid query length (exactly kmer_size)
        min_query = "QVQ"  # 3 characters, exactly kmer_size
        result = prefiltering.calculate_top_matches(min_query)
        assert result is not None
        assert result.query == min_query

        # Test longer query
        long_query = (
            "QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGRIIPILGIANYAQKFQGRVTITADKSTSTAYMELSSLRSEDTAVYYCAR"
        )
        result = prefiltering.calculate_top_matches(long_query)
        assert result is not None
        assert result.query == long_query

        shorter_than_kmer_size_query = "QV"
        with pytest.raises(ValueError):
            prefiltering.calculate_top_matches(shorter_than_kmer_size_query)


class TestPrefilteringEdgeCases:
    """Test edge cases and practical scenarios for both sequence types."""

    def test_small_gene_database_nucleotide(self):
        """Test behavior with small nucleotide gene database."""
        small_genes = {"SINGLE_GENE": "ACGTACGTACGTACGTACGTACGTACGTACGTACGT"}

        prefiltering = riot_na.Prefiltering(genes=small_genes, kmer_size=6, distance_threshold=1, top_n=5, modulo_n=1)

        result = prefiltering.calculate_top_matches("ACGTACGTACGTACGTACGTACGTACGTACGTACGT")
        assert result is not None
        assert len(result.top_matches) <= 1  # Can't have more matches than genes

    def test_small_gene_database_amino_acid(self):
        """Test behavior with small amino acid gene database."""
        small_genes = {
            "SINGLE_GENE": "QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGRIIPILGIANYAQKFQGRVTITADKSTSTAYMELSSLRSEDTAVYYCAR"
        }

        prefiltering = riot_na.Prefiltering(
            genes=small_genes, kmer_size=3, distance_threshold=1, top_n=5, modulo_n=1  # Appropriate for amino acids
        )

        result = prefiltering.calculate_top_matches(
            "QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGRIIPILGIANYAQKFQGRVTITADKSTSTAYMELSSLRSEDTAVYYCAR"
        )
        assert result is not None
        assert len(result.top_matches) <= 1  # Can't have more matches than genes

    def test_segment_vs_top_matches_consistency_nucleotide(self):
        """Test that segment matching gives same results as top_matches for single-domain nucleotide sequences."""
        genes = {
            "IGHV1-69*01": "CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTCACCAGCTATGCTATGCACTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAACCCTAACAGTGGTGGCACAAACTATGCACAGAAGTTTCAGGGCAGGGTCACCATGACCAGGGACACGTCCATCAGCACAGCCTACATGGAGCTGAGCAGGCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGA",
            "IGHV3-23*01": "GAAGTGCAGCTGCTGGAGTCTGGCGGAGGATTGGTACAGCCTGGCGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGCAGTTACTGGATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCAGCTATTAATGGTGGTAGCGGCACTTATTACGCTGACTCTGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTTTATCTGCAAATGAACAGCCTGAGAGCTGAGGACACGGCTGTGTATTATTGTGCTAAGGGTGAT",
            "IGKV1-27*01": "GACATCCAGATGACCCAGTCTCCATCCTCCCTGTCTGCATCTGTAGGAGACAGAGTCACCATCACTTGCCGGGCAAGTCAGATCAGTAGCAGCTATTTAAATTGGTATCAGCAGAAACCAGGGAAAGCCCCTAAGCTCCTGATCTATGCTGCATCCAGTTTGCAAAGTGGGGTCCCATCAAGGTTCAGTGGCAGTGGATCTGGGACAGATTTCACTCTCACCATCAGCAGTCTGCAACCTGAAGATTTTGCAACTTATTACTGTCAGCAATCTTATAGTACTCCTCTCACTTTCGGCGGAGGGACCAAGGTGGAGATCAAACGG",
        }

        prefiltering = riot_na.Prefiltering(
            genes=genes, kmer_size=6, distance_threshold=10, top_n=10, modulo_n=1, min_segment_length=20
        )

        # Test with a single-domain nucleotide sequence (should match IGHV1-69*01)
        single_domain_query = "CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTCACCAGCTATGCTATGCACTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAACCCTAACAGTGGTGGCACAAACTATGCACAGAAGTTTCAGGGCAGGGTCACCATGACCAGGGACACGTCCATCAGCACAGCCTACATGGAGCTGAGCAGGCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGA"

        # Get results from both methods
        top_matches_result = prefiltering.calculate_top_matches(single_domain_query)
        segment_result = prefiltering.calculate_segment_matches(single_domain_query)

        # Both should find results
        assert len(top_matches_result.top_matches) > 0, "Top matches should find at least one result"
        assert len(segment_result.segments) > 0, "Segment matching should find at least one segment"

        # For single domain, there should be exactly one segment that spans most/all of the query
        assert len(segment_result.segments) == 1, "Single domain should result in exactly one segment"
        segment = segment_result.segments[0]

        # The segment should span most of the query (allowing for some boundary variation)
        query_length = len(single_domain_query)
        segment_length = segment.query_end - segment.query_start
        assert (
            segment_length >= query_length * 0.8
        ), f"Segment should cover at least 80% of query length. Got {segment_length}/{query_length}"

        # The top gene matches should be consistent between methods
        # Extract gene IDs from top matches
        top_gene_ids = {match.gene_id for match in top_matches_result.top_matches}
        segment_gene_ids = {match.gene_id for match in segment.matching_genes}

        # There should be overlap in the top gene matches
        common_genes = top_gene_ids.intersection(segment_gene_ids)
        assert (
            len(common_genes) > 0
        ), f"Should have common genes between methods. Top: {top_gene_ids}, Segment: {segment_gene_ids}"

        # The best match should be the same (or very similar coverage)
        if top_matches_result.top_matches and segment.matching_genes:
            best_top_match = top_matches_result.top_matches[0]
            best_segment_match = max(segment.matching_genes, key=lambda x: x.coverage)

            # Best matches should have similar coverage (within reasonable tolerance)
            coverage_diff = abs(best_top_match.coverage - best_segment_match.coverage)
            max_coverage = max(best_top_match.coverage, best_segment_match.coverage)
            if max_coverage > 0:
                relative_diff = coverage_diff / max_coverage
                assert (
                    relative_diff <= 0.2
                ), f"Coverage should be similar. Top: {best_top_match.coverage}, Segment: {best_segment_match.coverage}"

    def test_segment_vs_top_matches_consistency_amino_acid(self):
        """Test that segment matching gives same results as top_matches for single-domain amino acid sequences."""
        genes = {
            "IGHV1-69*01": "QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGRIIPILGIANYAQKFQGRVTITADKSTSTAYMELSSLRSEDTAVYYCAR",
            "IGHV3-23*01": "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK",
            "IGKV1-27*01": "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPLT",
            "IGHV1-2*02": "QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGWINPNSGGTNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCAR",
            "IGKV1D-33*01": "DIQMTQSPSSLSASVGDRVTITCQASQDISNYLNWYQQKPGKAPKLLIYDASNLETGVPSRFSGSGSGTDFTFTISSLQPEDIATYYCQQYDNLP",
        }

        prefiltering = riot_na.Prefiltering(
            genes=genes,
            kmer_size=3,  # Appropriate for amino acids
            distance_threshold=4,
            top_n=10,
            modulo_n=1,
            min_segment_length=50,  # Shorter for amino acids
        )

        # Test with a single-domain amino acid sequence (should match IGHV1-2*02)
        single_domain_query = (
            "DIQMTQSPSSLSASVGDTVTITCQANGYLNWYQQRRGKAPKLLIYDGSKLERGVPSRFSGRRWGQEYNLTINNLQPEDIATYFCQVYEFVVPGTRLDLK"
        )

        # Get results from both methods
        top_matches_result = prefiltering.calculate_top_matches(single_domain_query)
        segment_result = prefiltering.calculate_segment_matches(single_domain_query)

        # Both should find results
        assert len(top_matches_result.top_matches) > 0, "Top matches should find at least one result"
        assert len(segment_result.segments) > 0, "Segment matching should find at least one segment"

        # For single domain, there should be exactly one segment that spans most/all of the query
        assert len(segment_result.segments) == 1, "Single domain should result in exactly one segment"
        segment = segment_result.segments[0]

        # The segment should span most of the query (allowing for some boundary variation)
        query_length = len(single_domain_query)
        segment_length = segment.query_end - segment.query_start
        assert (
            segment_length >= query_length * 0.8
        ), f"Segment should cover at least 80% of query length. Got {segment_length}/{query_length}"

        # The top gene matches should be consistent between methods
        # Extract gene IDs from top matches
        top_gene_ids = {match.gene_id for match in top_matches_result.top_matches[0:2]}
        segment_gene_ids = {match.gene_id for match in segment.matching_genes[0:2]}
        assert (
            top_gene_ids == segment_gene_ids
        ), f"Top gene ids should be the same. Top: {top_gene_ids}, Segment: {segment_gene_ids}"

        # There should be overlap in the top gene matches
        common_genes = top_gene_ids.intersection(segment_gene_ids)
        assert (
            len(common_genes) > 0
        ), f"Should have common genes between methods. Top: {top_gene_ids}, Segment: {segment_gene_ids}"

        # The best match should be the same (or very similar coverage)
        if top_matches_result.top_matches and segment.matching_genes:
            best_top_match = top_matches_result.top_matches[0]
            best_segment_match = max(segment.matching_genes, key=lambda x: x.coverage)

            assert (
                best_top_match.gene_id == best_segment_match.gene_id
            ), f"Best match should be the same. Top: {best_top_match.gene_id}, Segment: {best_segment_match.gene_id}"
            assert (
                best_top_match.coverage == best_segment_match.coverage
            ), f"Best match should have the same coverage. Top: {best_top_match.coverage}, Segment: {best_segment_match.coverage}"
