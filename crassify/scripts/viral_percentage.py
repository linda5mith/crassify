import pandas as pd
from tqdm import tqdm


class ViralPercentageCalculator:
    """
    Calculate % viral, % completeness, protein coverage, novelty score.
    """

    def __init__(self, matches, metadata, output_dir):
        self.matches = matches
        self.metadata = metadata
        self.output_dir = output_dir

    @staticmethod
    def most_frequent_virus(series):
        return series.mode().iloc[0] if not series.mode().empty else None

    def calculate_percentage_viral(self):
        pbar = tqdm(total=5, desc="Calculating % viral")

        # 1. Count how many proteins from each genome had a hit
        hits_per_genome = (
            self.matches.groupby("qseqid_genome_ID")["qseqid"]
            .nunique()
            .reset_index()
            .rename(columns={"qseqid": "protein_hits", "qseqid_genome_ID": "genome_ID"})
        )
        pbar.update(1)

        # 2. Pairwise aggregation (query vs subject)
        pairwise = (
            self.matches
            .groupby(["qseqid_genome_ID", "sseqid_genome_ID"])
            .agg({"length": "sum", "bitscore": "sum", 
                  "sseqid_genome_length": "first", 
                  "sseqid_virus": "first"})
            .reset_index()
            .rename(columns={
                "length": "total_aln_length",
                "bitscore": "total_bitscore",
                "sseqid_virus": "subject_virus"
            })
        )

        # 3. For each query genome, pick the subject with longest total protein alignment + highest bitscore
        best_hits = pairwise.loc[
            pairwise.groupby("qseqid_genome_ID")["total_bitscore"].idxmax(),
            ["qseqid_genome_ID", "sseqid_genome_ID", "subject_virus", "sseqid_genome_length", "total_aln_length"]
        ].rename(columns={
            "qseqid_genome_ID": "genome_ID",
            "sseqid_genome_ID": "top_species_hit_genome_accn",
            "subject_virus": "top_species_hit"
        })
        pbar.update(1)

        # 4. Merge metadata + counts + best hit info
        merged = (
            self.metadata
            .merge(hits_per_genome, on="genome_ID", how="left")
            .merge(best_hits, on="genome_ID", how="left")
        )
        merged.drop_duplicates(subset=["genome_ID"], inplace=True)

        # 5. Calculate viral/completeness/novelty stats
        merged["% contig viral"] = (
            (merged["total_aln_length"] * 3) / merged["genome_length"] * 100
        ).clip(upper=100)

        merged["% contig completeness"] = (
            merged["genome_length"] / merged["sseqid_genome_length"] * 100
        ).clip(upper=100)
        pbar.update(1)

        protein_counts = (
            self.metadata[self.metadata["input_seq"] == 1]
            .groupby("genome_ID")["protein_ID"]
            .nunique()
            .reset_index()
            .rename(columns={"protein_ID": "#_proteins"})
        )
        merged = merged.merge(protein_counts, on="genome_ID", how="left")

        merged["% proteins aligned"] = (
            merged["protein_hits"] / merged["#_proteins"] * 100
        )
        merged["novelty_score"] = (
            100 - merged["% proteins aligned"].fillna(0)
            - 0.5 * merged["% contig completeness"].fillna(0)
        )
        merged["is_novel"] = merged["novelty_score"] > 60
        pbar.update(2)

        print(merged.columns)

        # 6. Final output
        out = merged[[
            "genome_ID", "genome_length", "protein_hits", "top_species_hit","top_species_hit_genome_accn",
            "sseqid_genome_length", "total_aln_length", "% contig viral",
            "% contig completeness", "#_proteins", "% proteins aligned",
            "novelty_score", "is_novel"
        ]]
        out = out.round(3).sort_values(by=['protein_hits','% contig completeness'],ascending=False)
        out.to_csv(f"{self.output_dir}/percentage_viral.csv", index=False)
        pbar.close()
        return out

