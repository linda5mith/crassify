import pandas as pd
import numpy as np
import os
from tqdm import tqdm


class VirusDistanceCalculator:
    """
    Compute directional inter-genome distances and top hits.
    Keeps query (contig) → subject (reference) direction.
    """

    def __init__(self, matches_df, percentage_viral_df, out_dir):
        self.matches_df = matches_df
        self.percentage_viral_df = percentage_viral_df
        self.out_dir = out_dir

    def compute_distances(self):
        df = self.matches_df.copy()

        # Save any rows with missing genome IDs
        missing = df[(df.qseqid_genome_ID.isna()) | (df.sseqid_genome_ID.isna())]
        if not missing.empty:
            missing.to_csv(
                os.path.join(self.out_dir, "matches_missing_pairs.csv"), index=False
            )
            print(f"⚠️  Warning: {len(missing)} matches with missing genome IDs saved to matches_missing_pairs.csv")

        # Ensure IDs are strings
        df["qseqid_genome_ID"] = df["qseqid_genome_ID"].fillna("NA").astype(str)
        df["sseqid_genome_ID"] = df["sseqid_genome_ID"].fillna("NA").astype(str)

        # Aggregate per (query_genome → subject_genome)
        agg = df.groupby(["qseqid_genome_ID", "sseqid_genome_ID"]).agg({
            "conserved_AA_#": "sum",         # summed aligned AA
            "qseqid_genome_length": "first", # query genome length
            "sseqid_genome_length": "first", # subject genome length
            "length": "sum",                 # total alignment length
            "pident": "mean",                # avg % identity
            "sseqid_virus": "first"          # subject virus name
        }).reset_index()

        # Remove self-hits
        agg = agg[agg["qseqid_genome_ID"] != agg["sseqid_genome_ID"]]

        # Compute distance
        agg["avg_genome_length"] = (
            agg["qseqid_genome_length"] + agg["sseqid_genome_length"]
        ) / 2
        agg["distance"] = 1 - (agg["conserved_AA_#"] * 3) / agg["avg_genome_length"]

        # Round & reorder
        out = agg.sort_values(["qseqid_genome_ID", "distance"]).round(3)
        out = out[[
            "qseqid_genome_ID",
            "sseqid_genome_ID",
            "sseqid_virus",
            "distance",
            "length",   # → will rename
            "pident",   # → will rename
            "avg_genome_length",
            "qseqid_genome_length",
            "sseqid_genome_length"
        ]]
        out.rename(columns={
            "length": "total_aln_length",
            "pident": "avg_pid"
        }, inplace=True)

        # Save
        out.to_csv(os.path.join(self.out_dir, "distances.csv"), index=False)
        self.distances_df = out

        return out
