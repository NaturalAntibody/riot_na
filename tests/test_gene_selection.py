# pylint: disable=redefined-outer-name

import pytest

from riot_na.alignment.aa_gene_alignments import create_vjc_aligner_aa
from riot_na.alignment.gene_aligner import (
    create_aa_j_gene_aligner,
    create_aa_v_gene_aligner,
)
from riot_na.api import api_mp
from riot_na.api.riot_numbering import create_riot_aa
from riot_na.config import GENE_DB_DIR
from riot_na.data.model import InputType, Locus, Organism, Scheme

AA_GENES_DIR = GENE_DB_DIR / "gene_db" / "aa_genes_deduplicated"

HEAVY_QUERY = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDRGYTYGNFDHWGQGTLVTVSS"
KAPPA_QUERY = (
    "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPLTFGGGTKVEIK"
)

# Germline names that do not follow the plain IGHV<group>-<subgroup>*<allele> shape.
IRREGULAR_GENE_NAMES = [
    "IGHV1-69-2*01",
    "IGHV4-30-2*01",
    "IGHV7-4-1*01",
    "IGHV3-NL1*01",
    "IGHV3-30*i02",
    "IGKV1D-13*02",
    "IGHV2-70D*04",
    "IGLV2-14*01_g168t_t198g_c337t",
    "IGLV3-16*01_t145g_t146c",
]


@pytest.fixture(scope="module")
def riot_aa():
    return create_riot_aa(allowed_species=(Organism.HOMO_SAPIENS,))


@pytest.fixture(scope="module")
def vjc_aligner():
    return create_vjc_aligner_aa(allowed_species=(Organism.HOMO_SAPIENS,))


@pytest.fixture(scope="module")
def v_aligner():
    return create_aa_v_gene_aligner(allowed_species=(Organism.HOMO_SAPIENS,), aa_genes_dir=AA_GENES_DIR)


def test_align_to_selected_genes(vjc_aligner):
    # when
    alignments = vjc_aligner.produce_aa_alignments_for_genes(HEAVY_QUERY, v_gene="IGHV3-7*01", j_gene="IGHJ6*02")[0]

    # then
    assert alignments.v is not None
    assert alignments.v.target_id == "IGHV3-7*01"
    assert alignments.v.species == Organism.HOMO_SAPIENS
    assert alignments.v.locus == Locus.IGH

    assert alignments.j is not None
    assert alignments.j.target_id == "IGHJ6*02"

    # c genes are not assigned when germlines are selected by name
    assert alignments.c is None


def test_selected_genes_reproduce_searched_result(riot_aa):
    # given
    searched = riot_aa.run_on_sequence("header", HEAVY_QUERY, Scheme.IMGT)

    # when
    selected = riot_aa.run_on_sequence(
        "header", HEAVY_QUERY, Scheme.IMGT, v_gene=searched.v_call, j_gene=searched.j_call
    )

    # then
    assert selected == searched


def test_selected_genes_override_the_search(riot_aa):
    # when
    first = riot_aa.run_on_sequence("header", HEAVY_QUERY, Scheme.IMGT, v_gene="IGHV3-7*01", j_gene="IGHJ4*02")
    second = riot_aa.run_on_sequence("header", HEAVY_QUERY, Scheme.IMGT, v_gene="IGHV3-23*01", j_gene="IGHJ5*02")

    # then
    assert (first.v_call, first.j_call) == ("IGHV3-7*01", "IGHJ4*02")
    assert (second.v_call, second.j_call) == ("IGHV3-23*01", "IGHJ5*02")
    assert first.sequence_aa == second.sequence_aa == HEAVY_QUERY
    assert first.scheme_residue_mapping and second.scheme_residue_mapping

    # the germline changed, so the germline side of the alignment has to change too
    assert first.v_germline_alignment_aa != second.v_germline_alignment_aa
    assert first.j_germline_alignment_aa != second.j_germline_alignment_aa


@pytest.mark.parametrize("scheme", list(Scheme))
def test_selected_genes_number_in_every_scheme(riot_aa, scheme):
    # when
    result = riot_aa.run_on_sequence("header", HEAVY_QUERY, scheme, v_gene="IGHV3-23*01", j_gene="IGHJ4*02")

    # then
    assert result.v_call == "IGHV3-23*01"
    assert result.j_call == "IGHJ4*02"
    assert result.numbering_scheme == scheme.value
    assert result.scheme_residue_mapping
    assert result.exc is None


def test_light_chain_genes(riot_aa):
    # when
    result = riot_aa.run_on_sequence("header", KAPPA_QUERY, Scheme.IMGT, v_gene="IGKV1D-39*01", j_gene="IGKJ2*01")

    # then
    assert result.v_call == "IGKV1D-39*01"
    assert result.j_call == "IGKJ2*01"
    assert result.locus == Locus.IGK.value


@pytest.mark.parametrize("gene_name", IRREGULAR_GENE_NAMES)
def test_irregular_gene_names_are_resolved(v_aligner, gene_name):
    # when
    gene = v_aligner.get_gene(Organism.HOMO_SAPIENS, gene_name)

    # then
    assert gene.name == gene_name


def test_every_human_v_gene_is_selectable_by_name(v_aligner):
    """Every V germline in the database has to be reachable by its exact name."""
    # given
    genes = [gene for gene in v_aligner.gene_lookup.values() if gene.species == Organism.HOMO_SAPIENS]
    assert genes

    # then
    for gene in genes:
        assert v_aligner.get_gene(Organism.HOMO_SAPIENS, gene.name) is gene


@pytest.mark.parametrize("locus", list(Locus))
def test_every_human_j_gene_is_selectable_by_name(locus):
    """Every J germline in the database has to be reachable by its exact name."""
    # given
    j_aligner = create_aa_j_gene_aligner(organism=Organism.HOMO_SAPIENS, locus=locus, aa_genes_dir=AA_GENES_DIR)
    genes = list(j_aligner.gene_lookup.values())
    assert genes

    # then
    for gene in genes:
        assert j_aligner.get_gene(Organism.HOMO_SAPIENS, gene.name) is gene


def test_unknown_v_gene_is_rejected(riot_aa):
    with pytest.raises(AssertionError, match="IGHV3-7\\*99 not found"):
        riot_aa.run_on_sequence("header", HEAVY_QUERY, Scheme.IMGT, v_gene="IGHV3-7*99", j_gene="IGHJ4*02")


def test_unknown_j_gene_is_rejected(riot_aa):
    with pytest.raises(AssertionError, match="IGHJ9\\*01 not found"):
        riot_aa.run_on_sequence("header", HEAVY_QUERY, Scheme.IMGT, v_gene="IGHV3-23*01", j_gene="IGHJ9*01")


def test_j_gene_from_other_locus_is_rejected(riot_aa):
    # the J gene is looked up in the locus of the selected V gene
    with pytest.raises(AssertionError, match="IGKJ1\\*01 not found"):
        riot_aa.run_on_sequence("header", HEAVY_QUERY, Scheme.IMGT, v_gene="IGHV3-23*01", j_gene="IGKJ1*01")


def test_gene_of_other_species_is_rejected(riot_aa):
    with pytest.raises(AssertionError, match="not found in mouse"):
        riot_aa.run_on_sequence(
            "header",
            HEAVY_QUERY,
            Scheme.IMGT,
            v_gene="IGHV3-23*01",
            j_gene="IGHJ4*02",
            species=Organism.MUS_MUSCULUS,
        )


def test_both_genes_are_required(riot_aa):
    with pytest.raises(AssertionError, match="Both v_gene and j_gene"):
        riot_aa.run_on_sequence("header", HEAVY_QUERY, Scheme.IMGT, v_gene="IGHV3-23*01")


def test_gene_selection_rejected_in_multiple_domains_mode():
    # given
    riot = create_riot_aa(allowed_species=(Organism.HOMO_SAPIENS,), return_all_domains=True)

    # then
    with pytest.raises(AssertionError, match="multiple domains mode"):
        riot.run_on_sequence("header", HEAVY_QUERY, Scheme.IMGT, v_gene="IGHV3-23*01", j_gene="IGHJ4*02")


def test_file_run_validates_genes_before_starting_workers(tmp_path, monkeypatch):
    """
    Germline names have to be validated in the parent process. An exception raised in a
    multiprocessing pool initializer makes the pool respawn workers forever instead of failing.
    """

    # given
    input_fasta = tmp_path / "in.fasta"
    input_fasta.write_text(f">heavy\n{HEAVY_QUERY}\n")
    output_csv = tmp_path / "out.csv"

    class PoolStarted(Exception):
        pass

    def fail_on_pool(*_args, **_kwargs):
        raise PoolStarted

    monkeypatch.setattr(api_mp.mp, "Pool", fail_on_pool)

    def run(v_gene, j_gene):
        api_mp.run_on_file_mp(
            GENE_DB_DIR,
            input_fasta,
            output_csv,
            input_type=InputType.AA,
            allowed_species=(Organism.HOMO_SAPIENS,),
            v_gene=v_gene,
            j_gene=j_gene,
        )

    # then unknown germlines are rejected before any worker is started
    with pytest.raises(AssertionError, match="IGHJ9\\*99 not found"):
        run("IGHV3-23*01", "IGHJ9*99")

    # and known germlines pass validation, so the run gets as far as the worker pool
    with pytest.raises(PoolStarted):
        run("IGHV3-23*01", "IGHJ4*02")


def test_assert_genes_exist_does_not_align(vjc_aligner):
    # then
    vjc_aligner.assert_genes_exist("IGHV3-23*01", "IGHJ4*02")

    with pytest.raises(AssertionError, match="IGHV3-23\\*99 not found"):
        vjc_aligner.assert_genes_exist("IGHV3-23*99", "IGHJ4*02")
