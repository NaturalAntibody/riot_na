from pathlib import Path

import pandas as pd


def melt_heavy_light(df: pd.DataFrame) -> pd.DataFrame:
    df = df[["therapeutic", "heavy_sequence", "light_sequence"]].melt(
        ["therapeutic"], var_name="chain_type", value_name="sequence"
    )
    df["sequence_id"] = df["therapeutic"] + "_" + df["chain_type"].str.strip("_sequence")
    df = df.set_index("sequence_id")
    df = df[df["sequence"] != "na"]
    return df


DIR = Path(__file__).parent.parent.parent / "data" / "therapeutics"

if __name__ == "__main__":
    sequence_types = ["human", "humanized"]
    dataset = pd.read_parquet(DIR / "therapeutics.parquet")
    dataset = dataset.drop_duplicates(subset=["therapeutic"]).reset_index(drop=True)

    dataset_heavy_light = melt_heavy_light(dataset)

    dataset_heavy_light.to_csv(DIR / "therapeutics.csv")

    with open(DIR / "therapeutics_aa.fasta", "w") as file:
        for pdb_id_chain, sequence in dataset_heavy_light[["sequence"]].itertuples(name=None):
            file.write(f">{pdb_id_chain}\n")
            file.write(f"{sequence}\n")

    metadata = pd.read_csv(DIR / "metadata.csv", skiprows=1)
    metadata.columns = metadata.columns.str.strip()
    metadata = metadata[["therapeutic", "type"]].copy()
    metadata["therapeutic"] = metadata["therapeutic"].str.strip().str.lower().str.capitalize()
    metadata["type"] = metadata["type"].str.strip().str.lower()
    metadata_human = metadata[metadata["type"].isin(sequence_types)]

    dataset_human = dataset.merge(metadata_human[["therapeutic"]], on="therapeutic")

    dataset_heavy_light_human = melt_heavy_light(dataset_human)
    dataset_heavy_light_human.to_csv(DIR / "therapeutics_human.csv")

    with open(DIR / "therapeutics_human_aa.fasta", "w") as file:
        for pdb_id_chain, sequence in dataset_heavy_light_human[["sequence"]].itertuples(name=None):
            file.write(f">{pdb_id_chain}\n")
            file.write(f"{sequence}\n")
