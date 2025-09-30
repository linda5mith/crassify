from Bio import SeqIO
import pandas as pd
import os, re
from tqdm import tqdm


class GenomeProteomeMapper:
    """
    Map genome accessions to protein accessions and lengths from FASTA input.
    Produces a protein metadata DataFrame.
    """

    def __init__(self, path_to_prots, output_dir):
        self.path_to_prots = path_to_prots
        self.output_dir = output_dir

    def map_genome_to_proteome(self):
        """Handle directory of FASTAs or a single FASTA file."""
        if os.path.isdir(self.path_to_prots):
            prots = [
                os.path.join(root, f)
                for root, _, files in os.walk(self.path_to_prots)
                for f in files if f.endswith((".faa", ".fasta"))
            ]
            return self._process_protein_files(prots)
        elif os.path.isfile(self.path_to_prots):
            return self._process_single_file(self.path_to_prots)
        else:
            raise FileNotFoundError(f"No FASTA files found in {self.path_to_prots}")

    def _process_protein_files(self, prots):
        all_records = []
        for seq in tqdm(prots, desc="Mapping genome → protein IDs"):
            genome_ID = os.path.splitext(os.path.basename(seq))[0]
            for record in SeqIO.parse(seq, "fasta"):
                all_records.append({
                    "genome_ID": genome_ID,
                    "protein_ID": record.id,
                    "function": record.description.split(" ", 1)[-1],
                    "protein_length": len(record.seq)
                })
        df = pd.DataFrame(all_records)
        return self._add_genome_lengths(df)

    def _process_single_file(self, file_path):
        records = []
        for record in SeqIO.parse(file_path, "fasta"):
            header = record.description
            prot_len = len(record.seq)
            if "]" in header:
                match = re.search(r"\[(.*?)\]", header)
                genome = match.group(1) if match else None
                function = " ".join(header.split(" ")[1:]).split("[")[0].strip()
            else:
                genome = record.id.rsplit("_", 1)[0]
                function = None
            records.append({
                "genome_ID": genome,
                "protein_ID": record.id,
                "function": function,
                "protein_length": prot_len
            })
        df = pd.DataFrame(records)
        return self._add_genome_lengths(df)

    def _add_genome_lengths(self, df):
        """Approx genome length = sum protein lengths × 3."""
        genome_len_df = df.groupby("genome_ID")["protein_length"].sum().reset_index()
        genome_len_df["genome_length"] = genome_len_df["protein_length"] * 3
        df = df.merge(genome_len_df[["genome_ID", "genome_length"]], on="genome_ID", how="left")
        df["input_seq"] = 1
        return df
