import json
from typing import Optional

import pandas as pd

from riot_na.common.serialization_utils import base64_decode


def base64_decode_series(df, col):
    df = df[df[col].notna()].copy()
    df[col] = df[col].apply(base64_decode).apply(json.loads)
    return df


def not_na_df(df):
    return df[(df.ne("") & df.notna()).all(axis=1)]


def calculate_validation_flags_summary(df):
    df_notna = df[df["v_call"].notna() & (df["v_call"] != "") & df["productive"]]
    return df_notna["additional_validation_flags"].apply(pd.Series).apply(lambda x: x[x.notna()].value_counts()).T


def validation_flags_comparison(
    df_1,
    df_2,
    df_1_name: str = "old",
    df_2_name: str = "new",
    additional_properties: Optional[list[str]] = None,
):
    if additional_properties is None:
        additional_properties = ["productive", "complete_vdj"]
    df_1_notna = df_1[df_1["v_call"].notna() & (df_1["v_call"] != "") & df_1["productive"]]
    df_2_notna = df_2[df_2["v_call"].notna() & (df_2["v_call"] != "") & df_2["productive"]]
    common_index = df_1_notna.index.intersection(df_2_notna.index)
    flags_as_columns = df_1_notna.loc[common_index]["additional_validation_flags"].apply(pd.Series)

    for riot_property in additional_properties:
        flags_as_columns[riot_property] = df_1_notna.loc[common_index][riot_property]
    df_1_val_flags = flags_as_columns.apply(lambda x: x[x.notna()].value_counts()).T

    flags_as_columns = df_2_notna.loc[common_index]["additional_validation_flags"].apply(pd.Series)
    for riot_property in additional_properties:
        flags_as_columns[riot_property] = df_2_notna.loc[common_index][riot_property]
    df_2_val_flags = flags_as_columns.apply(lambda x: x[x.notna()].value_counts()).T

    comp = df_1_val_flags.join(df_2_val_flags, lsuffix=f"_{df_1_name}", rsuffix=f"_{df_2_name}")
    comp = comp.fillna(0)
    comp[f"False_{df_1_name}-{df_2_name}"] = comp[f"False_{df_1_name}"] - comp[f"False_{df_2_name}"]
    comp[f"True_{df_1_name}-{df_2_name}"] = comp[f"True_{df_1_name}"] - comp[f"True_{df_2_name}"]
    return comp[
        [
            f"False_{df_1_name}",
            f"False_{df_2_name}",
            f"False_{df_1_name}-{df_2_name}",
            f"True_{df_1_name}",
            f"True_{df_2_name}",
            f"True_{df_1_name}-{df_2_name}",
        ]
    ]
