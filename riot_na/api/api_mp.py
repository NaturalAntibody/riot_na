import itertools
import multiprocessing as mp
from pathlib import Path
from typing import Any, Optional

from Bio import SeqIO
from psutil import cpu_count
from tqdm import tqdm

from riot_na.api.riot_numbering import create_riot_aa, create_riot_nt
from riot_na.common.io import count_fasta_records, write_airr_iter_to_csv
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


class _WorkerNT:
    def __init__(
        self,
        allowed_species: Optional[list[Organism]] = None,
        scheme: Scheme = Scheme.IMGT,
        db_dir: Path = GENE_DB_DIR,
        return_all_domains: bool = False,
    ) -> None:
        self.numbering = create_riot_nt(
            allowed_species=allowed_species, db_dir=db_dir, return_all_domains=return_all_domains
        )
        self.scheme = scheme

    def __call__(self, fasta_record: SeqIO.SeqRecord) -> AirrRearrangementEntryNT | list[AirrRearrangementEntryNT]:
        try:
            res = self.numbering.run_on_sequence(
                fasta_record.description,
                str(fasta_record.seq),
                self.scheme,
            )
        except Exception as exc:
            print(fasta_record.description)
            print(self.scheme)
            print(f"seq: {fasta_record.seq}")
            raise exc

        return res


class _WorkerAA:
    def __init__(
        self,
        allowed_species: Optional[list[Organism]] = None,
        scheme: Scheme = Scheme.IMGT,
        db_dir: Path = GENE_DB_DIR,
        return_all_domains: bool = False,
        extend_alignment: bool = False,
    ) -> None:
        self.numbering = create_riot_aa(
            allowed_species=allowed_species, db_dir=db_dir, return_all_domains=return_all_domains
        )
        self.scheme = scheme
        self.extend_alignment = extend_alignment

    def __call__(self, fasta_record: SeqIO.SeqRecord) -> AirrRearrangementEntryAA | list[AirrRearrangementEntryAA]:
        try:
            res = self.numbering.run_on_sequence(
                fasta_record.description,
                str(fasta_record.seq),
                self.scheme,
                extend_alignment=self.extend_alignment,
            )
        except Exception as exc:
            print(fasta_record.description)
            print(self.scheme)
            print(f"seq: {fasta_record.seq}")
            raise exc

        return res


def _worker_initializer(
    allowed_species: Optional[list[Organism]] = None,
    scheme: Scheme = Scheme.IMGT,
    db_dir: Path = GENE_DB_DIR,
    input_type: InputType = InputType.NT,
    return_all_domains: bool = False,
    extend_alignment: bool = False,
):
    global worker  # pylint: disable=global-statement,global-variable-undefined
    if input_type == InputType.NT:
        worker = _WorkerNT(allowed_species, scheme, db_dir, return_all_domains)  # type: ignore
    else:
        worker = _WorkerAA(allowed_species, scheme, db_dir, return_all_domains, extend_alignment)  # type: ignore


def _worker_call(*args, **kwds):
    return worker(*args, **kwds)  # type: ignore


def run_on_file_mp(  # pylint: disable=too-many-arguments
    db_dir: Path,
    input_fasta_path: Path,
    result_path: Path,
    n_processes: int = cpu_count(logical=False),
    input_format: str = "fasta",
    scheme: Scheme = Scheme.IMGT,
    allowed_species: Optional[list[Organism]] = None,
    input_type: InputType = InputType.NT,
    limit: Optional[int] = None,
    return_all_domains: bool = False,
    extend_alignment: bool = False,
):
    with mp.Pool(
        processes=n_processes,
        initializer=_worker_initializer,
        initargs=(allowed_species, scheme, db_dir, input_type, return_all_domains, extend_alignment),
    ) as pool:
        result_iter = tqdm(
            pool.imap(_worker_call, itertools.islice(SeqIO.parse(input_fasta_path, input_format), limit)),
            total=count_fasta_records(input_fasta_path, input_format=input_format),
        )
        record_type: Any
        match input_type:
            case InputType.NT:
                record_type = AirrRearrangementEntryNT if not return_all_domains else SegmentedAirrRearrangementEntryNT
            case InputType.AA:
                record_type = AirrRearrangementEntryAA if not return_all_domains else SegmentedAirrRearrangementEntryAA

        write_airr_iter_to_csv(result_path, record_type, result_iter)


if __name__ == "__main__":
    # Example usage

    # .write.text("/home/pawel.dudzic/workspace/analyzer/projects/automation/experiments/exploration/new_therapeutics/paper/new_mess/data/13_08_2025/molecules_segments.fasta", lineSep="\n")
    INPUT_PATH = "/home/pawel.dudzic/workspace/analyzer/projects/automation/experiments/exploration/new_therapeutics/paper/new_mess/data/13_08_2025/molecules_segments.fasta"
    OUTPUT_PATH = "/home/pawel.dudzic/workspace/analyzer/projects/automation/experiments/exploration/new_therapeutics/paper/new_mess/data/13_08_2025/molecules_segments_numbered.csv"
    run_on_file_mp(
        db_dir=GENE_DB_DIR,
        input_fasta_path=Path(INPUT_PATH),
        result_path=Path(OUTPUT_PATH),
        n_processes=8,
        input_format="fasta",
        scheme=Scheme.IMGT,
        allowed_species=[Organism.HOMO_SAPIENS],
        input_type=InputType.AA,
    )
