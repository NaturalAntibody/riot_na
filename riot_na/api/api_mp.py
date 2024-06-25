import itertools
import multiprocessing as mp
from pathlib import Path
from typing import Optional

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
)


class _WorkerNT:
    def __init__(
        self,
        allowed_species: Optional[list[Organism]] = None,
        scheme: Scheme = Scheme.IMGT,
        db_dir: Path = GENE_DB_DIR,
    ) -> None:
        self.numbering = create_riot_nt(allowed_species=allowed_species, db_dir=db_dir)
        self.scheme = scheme

    def __call__(self, fasta_record: SeqIO.SeqRecord) -> AirrRearrangementEntryNT:
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
    ) -> None:
        self.numbering = create_riot_aa(allowed_species=allowed_species, db_dir=db_dir)
        self.scheme = scheme

    def __call__(self, fasta_record: SeqIO.SeqRecord) -> AirrRearrangementEntryAA:
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


def _worker_initializer(
    allowed_species: Optional[list[Organism]] = None,
    scheme: Scheme = Scheme.IMGT,
    db_dir: Path = GENE_DB_DIR,
    input_type: InputType = InputType.NT,
):
    global worker  # pylint: disable=global-statement,global-variable-undefined
    if input_type == InputType.NT:
        worker = _WorkerNT(allowed_species, scheme, db_dir)  # type: ignore
    else:
        worker = _WorkerAA(allowed_species, scheme, db_dir)  # type: ignore


def _worker_call(*args, **kwds):
    return worker(*args, **kwds)  # type: ignore


def run_on_file_mp(
    db_dir: Path,
    input_fasta_path: Path,
    result_path: Path,
    n_processes: int = cpu_count(logical=False),
    input_format: str = "fasta",
    scheme: Scheme = Scheme.IMGT,
    allowed_species: Optional[list[Organism]] = None,
    input_type: InputType = InputType.NT,
    limit: Optional[int] = None,
):
    with mp.Pool(
        processes=n_processes, initializer=_worker_initializer, initargs=(allowed_species, scheme, db_dir, input_type)
    ) as pool:
        result_iter = tqdm(
            pool.imap(_worker_call, itertools.islice(SeqIO.parse(input_fasta_path, input_format), limit)),
            total=count_fasta_records(input_fasta_path, input_format=input_format),
        )

        record_type = AirrRearrangementEntryNT if input_type == InputType.NT else AirrRearrangementEntryAA
        write_airr_iter_to_csv(result_path, record_type, result_iter)
