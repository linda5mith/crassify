import pandas as pd
import numpy as np
from tqdm import tqdm
import os


class VirusDistanceCalculator:
    """
    Compute inter-genome distances, top hits, and PHYLIP table.
    """

    def __init__(self, matches_df, percentage_viral_df, out_dir):
        self.matches_df = matches_df
        self.percentage_viral_df = percentage_viral_df
        self.out_dir = out_dir

    def compute_distances(self):
        df = self.matches_df.copy()

        df[(df.qseqid_genome_ID.isna()) | (df.sseqid_genome_ID.isna())].to_csv(
            os.path.join(self.out_dir, "matches_na_genomes.csv"), index=False
        )

        # Ensure IDs are strings and fill NaNs
        df["qseqid_genome_ID"] = df["qseqid_genome_ID"].fillna("NA").astype(str)
        df["sseqid_genome_ID"] = df["sseqid_genome_ID"].fillna("NA").astype(str)

        # Symmetric genome pairs
        def make_pair(a, b):
            try:
                return tuple(sorted([str(a), str(b)]))
            except Exception:
                return ("NA", "NA")

        df["genome_pair"] = df.apply(
            lambda r: make_pair(r["qseqid_genome_ID"], r["sseqid_genome_ID"]), axis=1
        )

        # Aggregate by pair
        agg = df.groupby("genome_pair").agg({
            "conserved_AA_#": "sum",
            "qseqid_genome_length": "first",
            "sseqid_genome_length": "first",
            "length": "sum",       # total alignment length
            "pident": "mean"       # average % identity across all alignments
        }).reset_index()

        # Identify bad pairs
        bad_pairs = agg[agg["genome_pair"].apply(lambda x: "NA" in x)]
        if not bad_pairs.empty:
            print(f"⚠️  Warning: {len(bad_pairs)} genome pairs contained NA values and were removed.")
            bad_pairs.to_csv(os.path.join(self.out_dir, "removed_na_pairs.csv"), index=False)

        # Remove self-hits and "NA" pairs
        agg = agg[agg["genome_pair"].apply(lambda x: x[0] != x[1])]
        agg = agg[~agg["genome_pair"].apply(lambda x: "NA" in x)]

        # Compute distances
        agg["avg_genome_length"] = (agg["qseqid_genome_length"] + agg["sseqid_genome_length"]) / 2
        agg["distance"] = 1 - (agg["conserved_AA_#"] * 3) / agg["avg_genome_length"]

        # Split genome_pair back
        agg[["qseqid_genome_ID", "sseqid_genome_ID"]] = pd.DataFrame(
            agg["genome_pair"].tolist(), index=agg.index
        )
        agg.drop(columns=["genome_pair"], inplace=True)

        # Save
        out = agg.sort_values(["qseqid_genome_ID", "distance"]).round(3)
        out = out[[
            "qseqid_genome_ID",
            "sseqid_genome_ID",
            "distance",
            "length",   # now = total_aln_length
            "pident",   # now = avg_pid
            "avg_genome_length",
            "qseqid_genome_length",
            "sseqid_genome_length"
        ]]
        out.rename(columns={"length": "total_aln_length", "pident": "avg_pid"}, inplace=True)

        out.to_csv(os.path.join(self.out_dir, "distances.csv"), index=False)

        self.distances_df = out
        return out


    # def create_phylip_table(self):
    #     if not hasattr(self, "distances_df"):
    #         raise ValueError("Run compute_distances first.")
    #     matrix = self.distances_df.pivot_table(
    #         index="qseqid_genome_ID", columns="sseqid_genome_ID", values="distance"
    #     )
    #     matrix = matrix.fillna(1).astype(float)
    #     np.fill_diagonal(matrix.values, 0)
    #     phylip_file = os.path.join(self.out_dir, "phylip.dist")
    #     matrix.to_csv(phylip_file, sep=" ")
    #     return phylip_file