# type: ignore
import json
import multiprocessing as mp
import re
from functools import partial

import pandas as pd
from tqdm import tqdm

from riot_na.alignment.alignment_utils import get_cigar_op_groups

GRID_SEARCH_PARAMS = ["top", "kmer_size", "distance_threshold", "modulo", "gap_open", "gap_ext", "x_drop"]


def parse_file_name(path, params=None):
    if params is None:
        params = GRID_SEARCH_PARAMS
    regex = "_".join([f"{param}_([0-9]*)" for param in params])
    match = re.search(regex, path.stem)
    return match.groups()


def _extract_gene_match(matches_string, top_n, with_rev_comp):
    matches = json.loads(matches_string)[:top_n]
    if with_rev_comp:
        return {(match["gene_id"], match["rev_comp"]) for match in matches}
    return {match["gene_id"] for match in matches}


def _agg_grid_search_worker(path, params, top_n, ground_truth_df, with_rev_comp):
    df = pd.read_csv(path, index_col=0, usecols=["sequence_id", "best_genes"])
    df["best_genes"] = df["best_genes"].apply(partial(_extract_gene_match, top_n=top_n, with_rev_comp=with_rev_comp))
    df = df.join(ground_truth_df, rsuffix="_true", how="inner")
    df["is_match"] = df.apply(lambda s: s["target_id"] in s["best_genes"], axis=1)
    match_percent = df["is_match"].value_counts(normalize=True).loc[True]
    return {**dict(zip(params, parse_file_name(path, params))), "match_percent": match_percent}


def aggregate_grid_search_results(results_dir, params, top_n, ground_truth_df, with_rev_comp):
    worker_partial = partial(
        _agg_grid_search_worker,
        params=params,
        top_n=top_n,
        ground_truth_df=ground_truth_df,
        with_rev_comp=with_rev_comp,
    )
    paths = list(results_dir.glob("top_*"))
    with mp.Pool(4) as pool:
        res = list(tqdm(pool.imap(worker_partial, paths), total=len(paths)))
        res_df = pd.DataFrame.from_records(res)
        res_df["top"] = top_n
        res_df[params] = res_df[params].astype(int)
        return res_df


def format_cigar(cigar: str, query: str, target: str):
    cigar_items = get_cigar_op_groups(cigar)
    out = ""
    pos = 0
    for cnt, operation in cigar_items:
        if operation == "M":
            for i in range(cnt):
                if query[pos + i] == target[pos + i]:
                    out += "|"
                else:
                    out += "X"
        elif operation == "=":
            out += cnt * "|"
        elif operation == "X":
            out += cnt * "X"
        elif operation == "I":
            target = target[:pos] + cnt * "-" + target[pos:]
            out += cnt * "+"
        elif operation == "D":
            query = query[:pos] + cnt * "-" + query[pos:]
            out += cnt * "-"
        else:
            out += cnt * "?"
        pos += cnt
    return out, query, target


def calculate_cigar_op_sum(cigar: str) -> int:
    return sum(len for len, _ in get_cigar_op_groups(cigar))


def display_riot_alignment(row):
    print(f"ID: {row.name}")
    print(f"rev_comp: {row['rev_comp']}")
    print()

    t_start = row["t_start"]
    t_end = row["t_end"]
    q_start = row["q_start"]
    q_end = row["q_end"]
    seq_identity = row["seq_identity"]
    bit_score = row["bit_score"]
    e_value = row["e_value"]
    target = row["target"][t_start:t_end]
    query = row["sequence"][q_start:q_end]
    cigar = row["cigar"]
    cigar_op_sum = calculate_cigar_op_sum(cigar)

    print(f"Riot: {row['target_id']}")
    print(f"{t_start=}, {t_end=}, {len(target)=}, {q_start=}, {q_end=}, {len(query)=}, {cigar=}, {cigar_op_sum=}")
    print(f"{seq_identity=}")
    print(f"{bit_score=}")
    print(f"{e_value=}")
    cigar_exploded, query, target = format_cigar(cigar, query, target)

    print("Target: ", target)
    print("        ", cigar_exploded)
    print("Query:  ", query)
