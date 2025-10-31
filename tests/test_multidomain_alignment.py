from typing import Dict, Optional

import pytest

from riot_na import Locus, Organism, Scheme
from riot_na.api.riot_numbering import (
    RiotNumberingAA,
    RiotNumberingNT,
    create_riot_aa,
    get_or_create_riot_aa,
    get_or_create_riot_nt,
)


class TestMultiDomainAlignments:
    """Test real-world multi-domain therapeutic antibody scenarios."""

    @pytest.fixture
    def human_riot_aa(self) -> RiotNumberingAA:
        """Get amino acid numbering API."""
        return create_riot_aa(allowed_species=[Organism.HOMO_SAPIENS], return_all_domains=True)

    @pytest.fixture
    def riot_aa(self) -> RiotNumberingAA:
        """Get amino acid numbering API."""
        return get_or_create_riot_aa(return_all_domains=True)

    @pytest.fixture
    def riot_nt(self) -> RiotNumberingNT:
        """Get nucleotide numbering API."""
        return get_or_create_riot_nt(return_all_domains=True)

    def test_devextinetug_amino_acid(self, riot_aa: RiotNumberingAA):
        """Test Devextinetug construct (amino acid)."""
        devextinetug_sequence = (
            "EIQLQQSGPELGKPGASVKVSCRASGFSFADYYIYWVKQSHGKSLELIGYIDPFNGGDTYNQIFKGKATLTVDKSSSTAFMYLNSLTSEDSAVYYCAAFRNPSFDFWGQGTTLTVSS"
            "GGGGGGGGSGGGGGGGGG"
            "QIVLIQSPPIMSASPGEKVTLTCSASSSVSSRYLYWYQQKPGSSPKLWIYGTSNLASGVPARFSGSGSGTSFSLTISSMEAEDAASYFCHQWSSFPFTFGSGTKLEIK"
        )

        expected_genes: Dict[str, Dict[str, Optional[str]]] = {
            "domain_1": {
                "v_call": "IGHV-3GZS",
                "j_call": "IGHJ-2FQV",
                "c_call": None,
                "locus": Locus.IGH.value,
            },
            "domain_2": {
                "v_call": "IGKV-3GON",
                "j_call": "IGKJ-Z5J4",
                "c_call": None,
                "locus": Locus.IGK.value,
            },
        }

        results = riot_aa.run_on_sequence("devextinetug", devextinetug_sequence, scheme=Scheme.IMGT)

        assert isinstance(results, list), "Should return list when return_all_domains=True"
        assert len(results) == len(expected_genes), f"Should detect {len(expected_genes)} domains"
        for i, domain in enumerate(results):
            assert (
                domain.locus == expected_genes[f"domain_{i+1}"]["locus"]
            ), f"Expected locus {expected_genes[f'domain_{i+1}']['locus']}, got {domain.locus}"

            assert (
                domain.v_call == expected_genes[f"domain_{i+1}"]["v_call"]
            ), f"Expected V call {expected_genes[f'domain_{i+1}']['v_call']}, got {domain.v_call}"
            assert (
                domain.j_call == expected_genes[f"domain_{i+1}"]["j_call"]
            ), f"Expected J call {expected_genes[f'domain_{i+1}']['j_call']}, got {domain.j_call}"

    def test_scfv_anti_cd19_amino_acid(self, riot_aa: RiotNumberingAA):
        """Test scFv construct (anti-CD19) - VH-linker-VL format (amino acid)."""
        # Real anti-CD19 scFv sequence with extended domains for better recognition
        # This represents a typical CAR-T or therapeutic scFv construct
        scfv_sequence = (
            # Extended VH domain (anti-CD19) - includes more of the constant region for better identification
            "EVQLVESGGGLVKPGGSLRLSCAASGFTFSNYGMHWVRQAPGKGLEWVGWINTYTGEPTYAADFKRRFTFSLDTSKSTAYLQMNSLRAEDTAVYYCAKYYGDYLDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"
            # (GGGGS)3 linker
            "GGGGSGGGGSGGGGS"
            # VL domain (anti-CD19)
            "EIVLTQSPGTLSLSPGERATLSCRASQRVSSSYLAWYQQKPGQAPRLLIYDASSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQQYGQSPYTFGQGTKVEIK"
        )

        expected_genes: Dict[str, Dict[str, Optional[str]]] = {
            "domain_1": {
                "v_call": "IGHV3-30*02",
                "j_call": "IGHJ4*02",
                "c_call": "IGHG1",
                "locus": Locus.IGH.value,
            },
            "domain_2": {
                "v_call": "IGKV3-20*01",
                "j_call": "IGKJ2*01",
                "c_call": None,
                "locus": Locus.IGK.value,
            },
        }

        # Test with return_all_domains=True to get all domains
        results = riot_aa.run_on_sequence("anti_CD19_scFv", scfv_sequence, scheme=Scheme.IMGT)

        # Should return multiple domains (expecting 2: VH and VL)
        assert isinstance(results, list), "Should return list when return_all_domains=True"
        assert len(results) == 2, "Should detect 2 domains"

        for i, domain in enumerate(results):
            assert (
                domain.locus == expected_genes[f"domain_{i+1}"]["locus"]
            ), f"Expected locus {expected_genes[f'domain_{i+1}']['locus']}, got {domain.locus}"
            assert (
                domain.v_call == expected_genes[f"domain_{i+1}"]["v_call"]
            ), f"Expected V call {expected_genes[f'domain_{i+1}']['v_call']}, got {domain.v_call}"
            assert (
                domain.j_call == expected_genes[f"domain_{i+1}"]["j_call"]
            ), f"Expected J call {expected_genes[f'domain_{i+1}']['j_call']}, got {domain.j_call}"
            assert (
                domain.c_call == expected_genes[f"domain_{i+1}"]["c_call"]
            ), f"Expected C call {expected_genes[f'domain_{i+1}']['c_call']}, got {domain.c_call}"

    def test_scfv_anti_her2_nucleotide(self, riot_nt: RiotNumberingNT):
        """Test scFv construct (anti-HER2) - nucleotide sequence."""
        # Real anti-HER2 scFv nucleotide sequence (trastuzumab-derived VH-linker-VL)
        scfv_nt_sequence = (
            # VH domain (nucleotide)
            "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTCCAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGCAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGCAAGGGGCTGGAGTGGGTGGCAGTTATATCATATGATGGAAGTAATAAATACTATGCAGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGAGATCGGGGGAGGAACGGTTATGATGCTTTTGATATCTGGGGCCAAGGGACAATGGTCACCGTCTCTTCAG"
            # (GGGGS)3 linker (nucleotide)
            "GGTGGCGGTTCAGGCGGAGGTGGCTCTGGCGGTGGCGGATCG"
            # VL domain (nucleotide)
            "GACATCCAGATGACCCAGTCTCCATCCTCCCTGTCTGCATCTGTAGGAGACAGAGTCACCATCACTTGCCGGGCAAGTCAGGGCATTAGCAGTTGGCTGGCCTGGTATCAGCAGAAACCAGGGAAAGCCCCTAAGCTCCTGATCTATGCTGCATCCAGTTTGCAAAGTGGGGTCCCATCAAGGTTCAGTGGCAGTGGATCTGGGACAGATTTCACTCTCACCATCAGCAGTCTGCAACCTGAAGATTTTGCAACTTACTACTGTCAGCAAGCTAACAGCTTCCCTTATACGTTCGGCCAAGGGACCAAGGTGGAAATCAAACGG"
        )

        expected_genes: Dict[str, Dict[str, Optional[str]]] = {
            "domain_1": {
                "v_call": "IGHV3-30*03",
                "d_call": "IGHD3-16*02",
                "j_call": "IGHJ3*02",
                "c_call": None,
                "locus": Locus.IGH.value,
            },
            "domain_2": {
                "v_call": "IGKV1-12*01",
                "d_call": None,
                "j_call": "IGKJ1*01",
                "c_call": None,
                "locus": Locus.IGK.value,
            },
        }

        # Test with different numbering schemes
        for scheme in [Scheme.IMGT, Scheme.CHOTHIA]:
            result = riot_nt.run_on_sequence("anti_HER2_scFv", scfv_nt_sequence, scheme=scheme)

            assert isinstance(result, list), "Should return list when return_all_domains=True"
            assert len(result) == 2, "Should detect 2 domains"

            for i, domain in enumerate(result):
                assert (
                    domain.locus == expected_genes[f"domain_{i+1}"]["locus"]
                ), f"Expected locus {expected_genes[f'domain_{i+1}']['locus']}, got {domain.locus}"
                assert (
                    domain.v_call == expected_genes[f"domain_{i+1}"]["v_call"]
                ), f"Expected V call {expected_genes[f'domain_{i+1}']['v_call']}, got {domain.v_call}"
                assert (
                    domain.d_call == expected_genes[f"domain_{i+1}"]["d_call"]
                ), f"Expected D call {expected_genes[f'domain_{i+1}']['d_call']}, got {domain.d_call}"
                assert (
                    domain.j_call == expected_genes[f"domain_{i+1}"]["j_call"]
                ), f"Expected J call {expected_genes[f'domain_{i+1}']['j_call']}, got {domain.j_call}"
                assert (
                    domain.c_call == expected_genes[f"domain_{i+1}"]["c_call"]
                ), f"Expected C call {expected_genes[f'domain_{i+1}']['c_call']}, got {domain.c_call}"

    def test_simple_scfv_with_constant_region(self, riot_aa: RiotNumberingAA):
        """Test a simpler scFv with part of constant region for better recognition."""
        # Single scFv with extended constant region (more realistic for API testing)
        scfv_with_constant = (
            # VH domain with partial constant region
            "QVQLQQWGAGLLKPSETLSLTCAVFGGSFSGYYWSWIRQPPGKGLEWIGEINHRGNTNDNPSLKSRVTISVDTSKNQFALKLSSVTAADTAVYYCARERGYTYGNFDHWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"
            # Linker to VL
            "GGGGSGGGGSGGGGS"
            # VL domain
            "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPLT"
        )

        expected_genes: Dict[str, Dict[str, Optional[str]]] = {
            "domain_1": {
                "v_call": "IGHV4-34*01",
                "j_call": "IGHJ1*01",
                "c_call": "IGHG1",
                "locus": Locus.IGH.value,
            },
            "domain_2": {
                "v_call": "IGKV1D-39*01",
                "j_call": None,
                "c_call": None,
                "locus": Locus.IGK.value,
            },
        }

        # Test with return_all_domains=True to detect multiple domains
        results = riot_aa.run_on_sequence("scFv_with_constant", scfv_with_constant, scheme=Scheme.IMGT)

        # Should return list with at least 1 domain
        assert isinstance(results, list), "Should return list when return_all_domains=True"
        assert len(results) == 2, "Should detect 2 domains"
        for i, domain in enumerate(results):
            assert (
                domain.locus == expected_genes[f"domain_{i+1}"]["locus"]
            ), f"Expected locus {expected_genes[f'domain_{i+1}']['locus']}, got {domain.locus}"
            assert (
                domain.v_call == expected_genes[f"domain_{i+1}"]["v_call"]
            ), f"Expected V call {expected_genes[f'domain_{i+1}']['v_call']}, got {domain.v_call}"
            assert (
                domain.j_call == expected_genes[f"domain_{i+1}"]["j_call"]
            ), f"Expected J call {expected_genes[f'domain_{i+1}']['j_call']}, got {domain.j_call}"
            assert (
                domain.c_call == expected_genes[f"domain_{i+1}"]["c_call"]
            ), f"Expected C call {expected_genes[f'domain_{i+1}']['c_call']}, got {domain.c_call}"

    def test_car_t_construct_with_multiple_domains(self, riot_aa: RiotNumberingAA):
        """Test CAR-T receptor construct with multiple functional domains."""
        # Realistic CAR-T construct: scFv + hinge + transmembrane + signaling domains
        # This tests very long multi-domain constructs used in cell therapy
        car_t_sequence = (
            # scFv domain (anti-tumor antigen)
            "QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGRIIPILGIANYAQKFQGRVTITADKSTSTAYMELSSLRSEDTAVYYCAR"
            "GGGGSGGGGSGGGGS"
            "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPLT"
            # CD8 hinge region
            "TTTPAPRPPTPAPTIASQPLSLRPEACRPAAGGAVHTRGLDFACD"
            # CD28 transmembrane domain
            "IYIWAPLAGTCGVLLLSLVITLYC"
            # CD28 costimulatory domain
            "KRGRKKLLYIFKQPFMRPVQTTQEEDGCSCRFPEEEEGGCEL"
            # CD3ζ signaling domain
            "RVKFSRSADAPAYQQGQNQLYNELNLGRREEYDVLDKRRGRDPEMGGKPRRKNPQEGLYNELQKDKMAEAYSEIGMKGERRRGKGHDGLYQGLSTATKDTYDALHMQALPPR"
        )

        expected_genes: Dict[str, Dict[str, Optional[str]]] = {
            "domain_1": {
                "v_call": "IGHV1-69*04",
                "j_call": None,
                "c_call": None,
                "locus": Locus.IGH.value,
            },
            "domain_2": {
                "v_call": "IGKV1D-39*01",
                "j_call": None,
                "c_call": None,
                "locus": Locus.IGK.value,
            },
        }

        # Test with IMGT scheme
        result = riot_aa.run_on_sequence("CAR_T_construct", car_t_sequence, scheme=Scheme.IMGT)

        assert isinstance(result, list), "Should return list when return_all_domains=True"
        assert len(result) == len(expected_genes), f"Should detect {len(expected_genes)} domains"
        for i, domain in enumerate(result):
            assert (
                domain.locus == expected_genes[f"domain_{i+1}"]["locus"]
            ), f"Expected locus {expected_genes[f'domain_{i+1}']['locus']}, got {domain.locus}"
            assert (
                domain.v_call == expected_genes[f"domain_{i+1}"]["v_call"]
            ), f"Expected V call {expected_genes[f'domain_{i+1}']['v_call']}, got {domain.v_call}"
            assert (
                domain.j_call == expected_genes[f"domain_{i+1}"]["j_call"]
            ), f"Expected J call {expected_genes[f'domain_{i+1}']['j_call']}, got {domain.j_call}"
            assert (
                domain.c_call == expected_genes[f"domain_{i+1}"]["c_call"]
            ), f"Expected C call {expected_genes[f'domain_{i+1}']['c_call']}, got {domain.c_call}"

    def test_diabody_construct(self, riot_aa: RiotNumberingAA):
        """Test diabody construct (VH-VL with short linker for dimerization)."""
        # Diabody format: VH1-short_linker-VL1 (forces heterodimerization)
        diabody_sequence = (
            # VH domain
            "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK"
            # Short linker (5 amino acids - forces inter-chain pairing)
            "GGGGGGGS"
            # VL domain
            "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPLT"
        )

        expected_genes: Dict[str, Dict[str, Optional[str]]] = {
            "domain_1": {
                "v_call": "IGHV3-23*01",
                "j_call": None,
                "c_call": None,
                "locus": Locus.IGH.value,
            },
            "domain_2": {
                "v_call": "IGKV1D-39*01",
                "j_call": None,
                "c_call": None,
                "locus": Locus.IGK.value,
            },
        }
        result = riot_aa.run_on_sequence("diabody_construct", diabody_sequence, scheme=Scheme.IMGT)

        # Should identify as heavy chain (VH first)
        assert isinstance(result, list), "Should return list when return_all_domains=True"
        assert len(result) == len(expected_genes), f"Should detect {len(expected_genes)} domains"
        for i, domain in enumerate(result):
            assert (
                domain.locus == expected_genes[f"domain_{i+1}"]["locus"]
            ), f"Expected locus {expected_genes[f'domain_{i+1}']['locus']}, got {domain.locus}"
            assert (
                domain.v_call == expected_genes[f"domain_{i+1}"]["v_call"]
            ), f"Expected V call {expected_genes[f'domain_{i+1}']['v_call']}, got {domain.v_call}"
            assert (
                domain.j_call == expected_genes[f"domain_{i+1}"]["j_call"]
            ), f"Expected J call {expected_genes[f'domain_{i+1}']['j_call']}, got {domain.j_call}"
            assert (
                domain.c_call == expected_genes[f"domain_{i+1}"]["c_call"]
            ), f"Expected C call {expected_genes[f'domain_{i+1}']['c_call']}, got {domain.c_call}"

    def test_trispecific_antibody_construct(self, riot_aa: RiotNumberingAA):
        """Test trispecific antibody construct with three different binding domains."""
        # Complex trispecific format: scFv1 + scFv2 + scFv3
        trispecific_sequence = (
            # First scFv (anti-target 1)
            "QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGRIIPILGIANYAQKFQGRVTITADKSTSTAYMELSSLRSEDTAVYYCAR"
            "GGGGSGGGGSGGGGS"
            "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPLT"
            # Long linker
            "GGGGSGGGGSGGGGSGGGGSGGGGS"
            # Second scFv (anti-target 2)
            "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK"
            "GGGGSGGGGSGGGGS"
            "EIVLTQSPGTLSLSPGERATLSCRASQRVSSSYLAWYQQKPGQAPRLLIYDASSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQQYGQSPYT"
            # Long linker
            "GGGGSGGGGSGGGGSGGGGSGGGGS"
            # Third scFv (anti-target 3)
            "QVQLQQWGAGLLKPSETLSLTCAVFGGSFSGYYWSWIRQPPGKGLEWIGEINHRGNTNDNPSLKSRVTISVDTSKNQFALKLSSVTAADTAVYYCAR"
            "GGGGSGGGGSGGGGS"
            "QSALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYDVSNRPSGVSNRFSGSKSGNTASLTISGLQAEDEADYYCSSYTSSSTLV"
        )

        expected_genes: Dict[str, Dict[str, Optional[str]]] = {
            "domain_1": {
                "v_call": "IGHV1-69*04",
                "j_call": None,
                "c_call": None,
                "locus": Locus.IGH.value,
            },
            "domain_2": {
                "v_call": "IGKV1D-39*01",
                "j_call": "IGKJ4*01",
                "c_call": None,
                "locus": Locus.IGK.value,
            },
            "domain_3": {
                "v_call": "IGHV3-23*01",
                "j_call": None,
                "c_call": None,
                "locus": Locus.IGH.value,
            },
            "domain_4": {
                "v_call": "IGKV3-20*01",
                "j_call": None,
                "c_call": None,
                "locus": Locus.IGK.value,
            },
            "domain_5": {
                "v_call": "IGHV4-34*01",
                "j_call": None,
                "c_call": None,
                "locus": Locus.IGH.value,
            },
            "domain_6": {
                "v_call": "IGLV2-14*03",
                "j_call": None,
                "c_call": None,
                "locus": Locus.IGL.value,
            },
        }

        result = riot_aa.run_on_sequence("trispecific_Ab", trispecific_sequence, scheme=Scheme.IMGT)

        assert isinstance(result, list), "Should return list when return_all_domains=True"
        assert len(result) == len(expected_genes), f"Should detect {len(expected_genes)} domains"
        for i, domain in enumerate(result):
            assert (
                domain.locus == expected_genes[f"domain_{i+1}"]["locus"]
            ), f"Expected locus for domain {i+1} {expected_genes[f'domain_{i+1}']['locus']}, got {domain.locus}"
            assert (
                domain.v_call == expected_genes[f"domain_{i+1}"]["v_call"]
            ), f"Expected V call for domain {i+1} {expected_genes[f'domain_{i+1}']['v_call']}, got {domain.v_call}"
            assert (
                domain.j_call == expected_genes[f"domain_{i+1}"]["j_call"]
            ), f"Expected J call for domain {i+1} {expected_genes[f'domain_{i+1}']['j_call']}, got {domain.j_call}"
            assert (
                domain.c_call == expected_genes[f"domain_{i+1}"]["c_call"]
            ), f"Expected C call for domain {i+1} {expected_genes[f'domain_{i+1}']['c_call']}, got {domain.c_call}"

    def test_trastuzumab_therapeutic_antibody(self, riot_aa: RiotNumberingAA):
        """Test a therapeutic antibody sequence (based on trastuzumab)."""
        # Realistic full-length antibody sequence (heavy chain with full constant region)
        # This represents what you'd see in actual therapeutic development
        trastuzumab_heavy = (
            # Variable region
            "EVQLVESGGGLVQPGGSLRLSCAASGFTFSNYGMHWVRQAPGKGLEWVGWINTYTGEPTYAADFKRRFTFSLDTSKSTAYLQMNSLRAEDTAVYYCAKYYGDYLDYWGQGTLVTVSS"
            # IgG1 constant region (partial for testing)
            "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKRVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"
        )

        expected_genes: Dict[str, Dict[str, Optional[str]]] = {
            "domain_1": {
                "v_call": "IGHV3-30*02",
                "j_call": "IGHJ4*02",
                "c_call": "IGHG1",
                "locus": Locus.IGH.value,
            },
        }

        result = riot_aa.run_on_sequence("trastuzumab_heavy", trastuzumab_heavy, scheme=Scheme.IMGT)

        assert isinstance(result, list), "Should return list when return_all_domains=True"
        assert len(result) == len(expected_genes), f"Should detect {len(expected_genes)} domains"
        for i, domain in enumerate(result):
            assert (
                domain.locus == expected_genes[f"domain_{i+1}"]["locus"]
            ), f"Expected locus for domain {i+1} {expected_genes[f'domain_{i+1}']['locus']}, got {domain.locus}"
            assert (
                domain.v_call == expected_genes[f"domain_{i+1}"]["v_call"]
            ), f"Expected V call for domain {i+1} {expected_genes[f'domain_{i+1}']['v_call']}, got {domain.v_call}"
            assert (
                domain.j_call == expected_genes[f"domain_{i+1}"]["j_call"]
            ), f"Expected J call for domain {i+1} {expected_genes[f'domain_{i+1}']['j_call']}, got {domain.j_call}"
            assert (
                domain.c_call == expected_genes[f"domain_{i+1}"]["c_call"]
            ), f"Expected C call for domain {i+1} {expected_genes[f'domain_{i+1}']['c_call']}, got {domain.c_call}"

    def test_acapatamab_therapeutic_antibody(self, riot_aa: RiotNumberingAA):
        """Test a therapeutic antibody sequence (based on acapatamab)."""

        acapatamab_full_sequence = "QVQLVESGGGLVKPGESLRLSCAASGFTFSDYYMYWVRQAPGKCLEWVAIISDGGYYTYYSDIIKGRFTISRDNAKNSLYLQMNSLKAEDTAVYYCARGFPLLRHGAMDYWGQGTLVTVSSGGGGSGGGGSGGGGSDIQMTQSPSSLSASVGDRVTITCKASQNVDTNVAWYQQKPGQAPKSLIYSASYVYWDVPSRFSGSASGTDFTLTISSVQSEDFATYYCQQYDQQLITFGCGTKLEIKSGGGGSEVQLVESGGGLVQPGGSLKLSCAASGFTFNKYAMNWVRQAPGKGLEWVARIRSKYNNYATYYADSVKDRFTISRDDSKNTAYLQMNNLKTEDTAVYYCVRHGNFGNSYISYWAYWGQGTLVTVSSGGGGSGGGGSGGGGSQTVVTQEPSLTVSPGGTVTLTCGSSTGAVTSGNYPNWVQQKPGQAPRGLIGGTKFLAPGTPARFSGSLLGGKAALTLSGVQPEDEAEYYCVLWYSNRWVFGGGTKLTVLGGGGDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPCEEQYGSTYRCVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGKGGGGSGGGGSGGGGSGGGGSGGGGSGGGGSDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPCEEQYGSTYRCVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"

        expected_genes: Dict[str, Dict[str, Optional[str]]] = {
            "domain_1": {
                "v_call": "IGHV-JRQS",
                "j_call": "IGHJ-EGQA",
                "c_call": None,
                "locus": Locus.IGH.value,
            },
            "domain_2": {
                "v_call": "IGKV1-16*02",
                "j_call": "IGKJ5*01",
                "c_call": None,
                "locus": Locus.IGK.value,
            },
            "domain_3": {
                "v_call": "IGHV-GNPX",
                "j_call": "IGHJ-KZ3G",
                "c_call": None,
                "locus": Locus.IGH.value,
            },
            "domain_4": {
                "v_call": "IGLV7-43*01",
                "j_call": "IGLJ3*02",
                "c_call": "IGLC6",
                "locus": Locus.IGL.value,
            },
        }

        result = riot_aa.run_on_sequence("acapatamab_heavy", acapatamab_full_sequence, scheme=Scheme.IMGT)

        assert isinstance(result, list), "Should return list when return_all_domains=True"
        assert len(result) == len(expected_genes), f"Should detect {len(expected_genes)} domains"
        for i, domain in enumerate(result):
            domain_key = f"domain_{i+1}"
            expected = expected_genes[domain_key]
            assert (
                domain.locus == expected["locus"]
            ), f"Expected locus for {domain_key} {expected['locus']}, got {domain.locus}"
            assert (
                domain.v_call == expected["v_call"]
            ), f"Expected V call for {domain_key} {expected['v_call']}, got {domain.v_call}"
            assert (
                domain.j_call == expected["j_call"]
            ), f"Expected J call for {domain_key} {expected['j_call']}, got {domain.j_call}"
            assert (
                domain.c_call == expected["c_call"]
            ), f"Expected C call for {domain_key} {expected['c_call']}, got {domain.c_call}"
