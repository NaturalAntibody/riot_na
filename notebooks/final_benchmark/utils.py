from typing import Literal

import pandas as pd


def calculate_gene_allele_assignment_precision(
    ground_truth_df: pd.DataFrame, df: pd.DataFrame, gene: Literal["v", "j", "d", "c"] = "v"
) -> dict[str, float]:
    joined = ground_truth_df[[f"{gene}_call"]].join(df[[f"{gene}_call"]], lsuffix="_true", how="inner")
    return {
        "gene": (joined[f"{gene}_call_true"].str.split("*").str[0] == joined[f"{gene}_call"].str.split("*").str[0])
        .value_counts(normalize=True)
        .loc[True]
        * 100,
        "allele": (joined[f"{gene}_call_true"] == joined[f"{gene}_call"]).value_counts(normalize=True).loc[True] * 100,
    }
