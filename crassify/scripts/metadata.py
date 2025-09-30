import pandas as pd
from tqdm import tqdm
import re


def strip_version(pid):
    """Remove version suffix from protein IDs (e.g. QWM91148.2 -> QWM91148)."""
    return re.sub(r"\.\d+$", "", str(pid))


class MetadataMerger:
    """
    Merge generated protein metadata with DIAMOND matches.
    Handles protein_IDs with or without version suffix.
    """

    def __init__(self, matches, generated_metadata, output_dir):
        self.matches = matches
        self.generated_metadata = generated_metadata
        self.output_dir = output_dir
        self.metadata = None

    def load_and_concat_metadata(self, input_metadata_path):
        input_metadata = pd.read_csv(input_metadata_path, low_memory=False)

        # concat with generated metadata
        self.metadata = pd.concat(
            [self.generated_metadata, input_metadata], ignore_index=True
        )

        # normalize protein IDs (strip version)
        self.metadata["protein_ID_base"] = self.metadata["protein_ID"].map(strip_version)
        self.metadata.drop_duplicates(subset=["protein_ID_base"], inplace=True)

        # fill missing sequence flag
        self.metadata["input_seq"] = self.metadata["input_seq"].fillna(0)

    def merge_metadata_matches(self):
        matches = pd.read_csv(self.matches, sep="\t", header=None, low_memory=False)
        matches.columns = [
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"
        ]
        matches = matches[matches["qseqid"] != matches["sseqid"]]

        # normalize IDs in matches too
        matches["qseqid_base"] = matches["qseqid"].map(strip_version)
        matches["sseqid_base"] = matches["sseqid"].map(strip_version)

        return self._merge_with_metadata(matches)

    def _merge_with_metadata(self, matches):
        pbar = tqdm(total=2, desc="Merging metadata")

        # qseqid merge
        matches = matches.merge(
            self.metadata,
            left_on="qseqid_base",
            right_on="protein_ID_base",
            how="left"
        )
        matches.rename(columns={
            "genome_length": "qseqid_genome_length",
            "genome_ID": "qseqid_genome_ID",
            "function": "qseqid_function",
            "protein_length": "qseqid_protein_length",
            "virus": "qseqid_virus",
            "input_seq": "qseqid_sequence"
        }, inplace=True)
        matches.drop(columns=["protein_ID", "protein_ID_base"], inplace=True, errors="ignore")
        pbar.update(1)

        # sseqid merge
        matches = matches.merge(
            self.metadata,
            left_on="sseqid_base",
            right_on="protein_ID_base",
            how="left",
            suffixes=("", "_s")
        )
        matches.rename(columns={
            "genome_length": "sseqid_genome_length",
            "genome_ID": "sseqid_genome_ID",
            "function": "sseqid_function",
            "protein_length": "sseqid_protein_length",
            "virus": "sseqid_virus",
            "input_seq": "sseqid_sequence"
        }, inplace=True)
        matches.drop(columns=["protein_ID", "protein_ID_base", "input_seq_s", "virus_s"], inplace=True, errors="ignore")
        pbar.update(1)

        # conserved AA count
        matches["conserved_AA_#"] = (matches["length"] * matches["pident"]) / 100
        pbar.close()

        matches = matches.dropna(axis=1, how="all")
        return matches
