import time
from pathlib import Path
from typing import Any, Optional

import click
import psutil

from riot_na.api.api_mp import run_on_file_mp
from riot_na.api.riot_numbering import create_riot_aa, create_riot_nt
from riot_na.common.io import write_airr_iter_to_csv
from riot_na.config import GENE_DB_DIR
from riot_na.data.model import (
    AirrRearrangementEntryAA,
    AirrRearrangementEntryNT,
    InputType,
    Organism,
    Scheme,
    SegmentedAirrRearrangementEntryAA,
    SegmentedAirrRearrangementEntryNT,
)


@click.command()
@click.option("-f", "--input-file", type=Path, help="Path to input FASTA file.")
@click.option("-s", "--sequence", type=str, help="Input sequence file.")
@click.option("-o", "--output-file", type=Path, help="Path to output CSV file. If not specified, stdout is used.")
@click.option(
    "--scheme",
    type=click.Choice(Scheme),  # type: ignore
    default=Scheme.IMGT,
    help="Which numbering scheme should be used: imgt, kabat, chothia, martin. Default IMGT",
)
@click.option(
    "--species",
    type=click.Choice(Organism),  # type: ignore
    default=None,
    help="Which species germline sequences should be used. Default is all species.",
)
@click.option(
    "--input-type",
    type=click.Choice(InputType),  # type: ignore
    default=InputType.NT,
    help="What kind of sequences are provided on input. Default is nucleotide sequences.",
)
@click.option(
    "-p",
    "--ncpu",
    type=int,
    default=psutil.cpu_count(logical=False),
    help="Number of parallel processes to use. Default is number of physical cores.",
)
@click.option(
    "-e",
    "--extend_alignment",
    type=bool,
    default=False,
    help=(
        "Include unaligned beginning of the query sequence in numbering."
        "This option impacts only amino acid sequences."
    ),
)
@click.option(
    "--multiple-domains",
    type=bool,
    default=False,
    help=("Return all domains of multiple domain proteins."),
)
def run_riot(
    input_file: Optional[Path],
    sequence: Optional[str],
    output_file: Optional[Path],
    scheme: Scheme,
    species: Optional[Organism],
    input_type: InputType,
    ncpu: int,
    extend_alignment: bool,
    multiple_domains: bool,
):
    species_list = [species] if species else None
    if input_file and input_file.exists():
        if not output_file:
            output_file = Path() / "numbering_result.csv"
        start = time.perf_counter()
        run_on_file_mp(
            GENE_DB_DIR,
            input_file,
            output_file,
            input_format="fasta",
            scheme=scheme,
            allowed_species=species_list,
            n_processes=ncpu,
            input_type=input_type,
            return_all_domains=multiple_domains,
        )
        end = time.perf_counter()
        elapsed_time = end - start
        print("Execution time:", time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    elif sequence:
        result: Any
        record_type: Any
        match input_type:
            case InputType.NT:
                numbering_nt = create_riot_nt(
                    allowed_species=species_list, db_dir=GENE_DB_DIR, return_all_domains=multiple_domains
                )
                record_type = AirrRearrangementEntryNT if not multiple_domains else SegmentedAirrRearrangementEntryNT
                result = numbering_nt.run_on_sequence(header="-", query_sequence=sequence, scheme=scheme)
            case InputType.AA:
                numbering_aa = create_riot_aa(
                    allowed_species=species_list, db_dir=GENE_DB_DIR, return_all_domains=multiple_domains
                )
                record_type = AirrRearrangementEntryAA if not multiple_domains else SegmentedAirrRearrangementEntryAA
                result = numbering_aa.run_on_sequence(
                    header="-", query_sequence=sequence, scheme=scheme, extend_alignment=extend_alignment
                )

        if output_file:
            write_airr_iter_to_csv(output_file, record_type, result if not multiple_domains else [result])
        else:
            if not multiple_domains:
                print(result.__dict__)
            else:
                for single_result in result:
                    print(single_result.__dict__)

    else:
        print("Need to specify input sequence or FASTA file!")
