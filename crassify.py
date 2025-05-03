from Bio import SeqIO
import pandas as pd
import numpy as np
import os
import re
import argparse
import yaml
import subprocess
from tqdm import tqdm
from rich.console import Console
from rich.markdown import Markdown
from datetime import datetime
import difflib

class CircularityChecker:
    @staticmethod
    def find_circular_overlap(seq, min_overlap=100, max_overlap=1000, match_threshold=0.9):
        """
        Detects potential circularity in a nucleotide sequence by checking
        similarity between the start and end regions.
        """
        for overlap in range(max_overlap, min_overlap - 1, -1):
            start = seq[:overlap]
            end = seq[-overlap:]
            similarity = difflib.SequenceMatcher(None, start, end).ratio()
            if similarity >= match_threshold:
                return overlap
        return None

class GenomeProteomeMapper:
    def __init__(self, path_to_prots, output_dir):
        self.path_to_prots = path_to_prots
        self.output_dir = output_dir

    def map_genome_to_proteome(self):
        '''Map genome accession to protein accession and return genome_length, protein_length'''
        if os.path.isdir(self.path_to_prots):
            prots = [os.path.join(root, f) for root, _, files in os.walk(self.path_to_prots) for f in files if 'gff' not in f]
            if prots:
                return self._process_protein_files(prots)
        elif os.path.isfile(self.path_to_prots):
            return self._process_single_file(self.path_to_prots)
        else:
            print(f'{self.path_to_prots} is empty or does not contain fasta files.\nCheck -p path/to/proteins.')

    def _process_protein_files(self, prots):
        '''Takes a folder of protein fasta files as input, maps the genome to protein ID's and gets the approximate genome length.'''
        try:
            genome_data = {}
            for seq in tqdm(prots, desc='Mapping genome ID to protein IDs'):
                genome_ID = os.path.splitext(os.path.basename(seq))[0]
                records = list(SeqIO.parse(seq, 'fasta'))

                record_data = [{
                    'protein_ID': record.id,
                    'protein_length': len(record.seq)
                } for record in records]

                genome_len = sum(len(record.seq) for record in records) * 3

                circular_overlap = None
                for record in records:
                    circular_overlap = CircularityChecker.find_circular_overlap(str(record.seq))
                    if circular_overlap:
                        print(f'{record} IS CIRCULAR!!')
                        break
                genome_data[genome_ID] = {
                    'protein_data': record_data,
                    'genome_length': genome_len,
                    'is_circular': bool(circular_overlap),
                    'circular_overlap_length': circular_overlap
                }

            return self._save_metadata(genome_data)
        except Exception as e:
            print(e)

    def _process_single_file(self, file_path):
        try:
            genome, protein_IDs, protein_func, protein_lengths = [], [], [], []
            total_records = sum(1 for _ in SeqIO.parse(file_path, 'fasta'))
            with open(file_path, 'r') as file:
                for record in tqdm(SeqIO.parse(file, 'fasta'), total=total_records, desc='Mapping genome ID to protein ID\'s'):
                    self._process_record(record, genome, protein_IDs, protein_func, protein_lengths)
            
            df = self._create_dataframe(genome, protein_IDs, protein_func, protein_lengths)
            df.to_csv(os.path.join(self.output_dir, 'protein_metadata.csv'), index=False)
            empress_metadata = df
            empress_metadata['genome_ID'] = empress_metadata['genome_ID'].str.replace(' ', '_')
            empress_metadata.to_csv(os.path.join(self.output_dir, 'feature_metadata.tsv'), index=False, sep='\t')
            return df
        except Exception as e:
            print(e)
    
    def _process_record(self, record, genome, protein_IDs, protein_func, protein_lengths):
        protein_ID = record.id
        header = record.description
        prot_len = len(record.seq)
        if ']' in header:
            definition = re.search(r'\[(.*?)\]', header)
            if definition:
                definition = definition.group(1)
            else:
                definition = None
            function = ' '.join(header.split(' ')[1:]).split('[')[0].strip()
            genome.append(definition)
            protein_IDs.append(protein_ID)
            protein_func.append(function)
            protein_lengths.append(prot_len)  # Store protein lengths
        else:
            contig = protein_ID.rsplit('_', 1)[0]
            genome.append(contig)
            protein_func.append(None)
            protein_IDs.append(protein_ID)
            protein_lengths.append(prot_len)  # Store protein lengths

    def _create_dataframe(self, genome, protein_IDs, protein_func, protein_lengths):
        df = pd.DataFrame({
            'genome_ID': genome,
            'protein_ID': protein_IDs,
            'function': protein_func,
            'protein_length': protein_lengths  # Include protein lengths in the DataFrame
        })
        
        # Calculate genome lengths
        genome_length_df = df.groupby('genome_ID')['protein_length'].sum().reset_index()
        genome_length_df['genome_length'] = genome_length_df['protein_length'] * 3  # Assuming each protein length is in base pairs
        
        # Merge genome length with the main DataFrame
        df = df.merge(genome_length_df[['genome_ID', 'genome_length']], on='genome_ID', how='left')
        
        df['input_seq'] = 1
        return df


class MetadataMerger:
    def __init__(self, matches, generated_metadata, output_dir):
        self.matches = matches
        self.generated_metadata = generated_metadata  # Metadata generated by the program
        self.output_dir = output_dir
        self.metadata = None  # This will hold the concatenated metadata

    def load_and_concat_metadata(self, input_metadata_path):
        """
        Loads metadata from the input file and concatenates it with the metadata
        generated by the genome_proteome_mapper.
        """
        input_metadata = pd.read_csv(input_metadata_path, low_memory=False)

        # Concatenate the input DB metadata and generated metadata
        self.metadata = pd.concat([self.generated_metadata, input_metadata], ignore_index=True)
        self.metadata = self.metadata.drop_duplicates(subset=['protein_ID'], keep='first')
        self.metadata['input_seq'] = self.metadata['input_seq'].fillna(0)

    def merge_metadata_matches(self):
        '''Merges protein metadata with DIAMOND output matches.tsv'''
        matches = pd.read_csv(self.matches, sep='\t', header=None, low_memory=False)
        if matches.empty:
            print("No hits found between inputted proteins. Exiting program.")
            sys.exit(1)

        matches.columns = ['qseqid', 'sseqid', "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
        matches = matches[matches['qseqid'] != matches['sseqid']]
        matches = self._merge_with_metadata(matches)
        matches.to_csv(os.path.join(self.output_dir, 'matches.csv'), index=False)
        return matches

    def _merge_with_metadata(self, matches):
        if self.metadata is None:
            raise ValueError("Metadata has not been loaded or concatenated. Call load_and_concat_metadata first.")
        total_steps = 4
        pbar = tqdm(total=total_steps, desc='Merging metadata...')
        # Merge for qseqid
        matches = matches.merge(
            self.metadata, 
            left_on=['qseqid'], 
            right_on=['protein_ID'], 
            how='left', 
            suffixes=('', '_qseqid')
        ).rename(columns={
            'genome_length': 'qseqid_genome_length',
            'protein_ID': 'qseqid_ID',
            'genome_ID': 'qseqid_genome_ID',
            'function': 'qseqid_function',
            'virus': 'qseqid_virus',
            'input_seq': 'qseqid_input_seq',
            'protein_length': 'qseqid_protein_length'
        })
        pbar.update(1)
        # Merge for sseqid
        matches = matches.merge(
            self.metadata, 
            left_on=['sseqid'], 
            right_on=['protein_ID'], 
            how='left', 
            suffixes=('', '_sseqid')
        ).rename(columns={
            'genome_length': 'sseqid_genome_length',
            'protein_ID': 'sseqid_ID',
            'genome_ID': 'sseqid_genome_ID',
            'function': 'sseqid_function',
            'virus': 'sseqid_virus',
            'input_seq': 'sseqid_input_seq',
            'protein_length': 'sseqid_protein_length'
        })
        pbar.update(1)
        # Fill NaNs for input sequences
        matches[['qseqid_input_seq', 'sseqid_input_seq']] = matches[['qseqid_input_seq', 'sseqid_input_seq']].fillna(0)
        pbar.update(1)
        # Calculate conserved amino acids
        matches['conserved_AA_#'] = (matches['length'] / 100) * matches['pident'].abs()
        pbar.update(1)
        return matches


class ViralPercentageCalculator:
    def __init__(self, matches, metadata, output_dir):
        self.matches = matches
        self.metadata = metadata
        self.output_dir = output_dir

        self.matches.to_csv(os.path.join(self.output_dir, 'matches_%_calculator.csv'))
        self.metadata.to_csv(os.path.join(self.output_dir, 'metadata_%_calculator.csv'))

    @staticmethod
    def most_frequent_virus(series):
        return series.mode().iloc[0] if not series.mode().empty else None
    
    def calculate_percentage_viral(self):
        total_steps = 5
        pbar = tqdm(total=total_steps, desc='Calculating % viral for input contigs...')
        
        # 1: Calculate hits per genome
        hits_per_genome = self.matches.groupby('qseqid_genome_ID')['qseqid_ID'].nunique().reset_index()
        hits_per_genome.columns = ['genome_ID', 'protein_hits']
        pbar.update(1)

        # 2: Aggregate protein lengths and find the most frequent sseqid_virus
        aggregated = self.matches.groupby('qseqid_genome_ID').agg({
            'qseqid_protein_length': 'sum',
            'sseqid_protein_length': 'sum',
            'sseqid_genome_length': 'first',
            'sseqid_virus': self.most_frequent_virus
        }).reset_index()
        pbar.update(1)
        
        # 3: Merge hits and aggregates with metadata
        merged = pd.merge(self.metadata, hits_per_genome, on='genome_ID', how='left')
        merged = pd.merge(merged, aggregated, left_on='genome_ID', right_on='qseqid_genome_ID', how='left')
        merged = merged.drop_duplicates(subset=['genome_ID'], keep='first')
        merged['sseqid_protein_length'] = pd.to_numeric(merged['sseqid_protein_length'], errors='coerce')
        merged['genome_length'] = pd.to_numeric(merged['genome_length'], errors='coerce')

        # 4: Calculate % viral and % completeness
        merged['% contig viral'] = ((merged['sseqid_protein_length'] * 3) / merged['genome_length'] * 100).round(2).clip(upper=100)
        merged['% contig completeness'] = (merged['genome_length'] / merged['sseqid_genome_length'] * 100).round(2).clip(upper=100)
        pbar.update(1)

        # 5: Calculate protein counts and novelty score
        protein_counts = self.metadata[self.metadata['input_seq'] == 1].groupby('genome_ID')['protein_ID'].nunique().reset_index()
        protein_counts.columns = ['genome_ID', '#_proteins']
        merged = pd.merge(merged, protein_counts, on='genome_ID', how='left')

        merged['% proteins aligned'] = (merged['protein_hits'] / merged['#_proteins'] * 100).round(2)
        merged['novelty_score'] = (
            100 
            - merged['% proteins aligned'].fillna(0)
            - 0.5 * merged['% contig completeness'].fillna(0)
        ).round(2)
        merged['is_novel'] = merged['novelty_score'] > 60
        pbar.update(1)

        # Final output
        merged = merged[['genome_ID', 'genome_length', 'protein_hits', 'sseqid_virus', 'sseqid_genome_length',
                        'sseqid_protein_length', '% contig viral', '% contig completeness', '#_proteins',
                        '% proteins aligned', 'novelty_score', 'is_novel']]

        merged = merged.rename(columns={
            'sseqid_protein_length': 'total_aln_length',
            'sseqid_virus': 'top_species_hit'
        })

        merged = merged.sort_values(by=['novelty_score', 'genome_length'], ascending=[False, False])
        merged.to_csv(os.path.join(self.output_dir, 'percentage_viral.csv'), index=False)
        pbar.close()
        return merged

class VirusDistanceCalculator:
    def __init__(self, matches_df, percentage_viral_df, out_dir):
        """
        Initialize the VirusDistanceCalculator with matches and percentage_viral DataFrames.
        
        :param matches_df: DataFrame containing match data.
        :param percentage_viral_df: DataFrame containing percentage viral data.
        """
        self.matches_df = matches_df
        self.percentage_viral_df = percentage_viral_df
        self.out_dir = out_dir

    def filter_viruses(self, viral_threshold=50, completeness_threshold=70):
        """
        Filter viruses based on viral and completeness thresholds.
        
        :param viral_threshold: Minimum % contig viral to include.
        :param completeness_threshold: Minimum % contig completeness to include.
        :return: Filtered DataFrame.
        """
        self.filtered_viruses = self.percentage_viral_df[
            (self.percentage_viral_df['% contig viral'] > viral_threshold) & 
            (self.percentage_viral_df['% contig completeness'] > completeness_threshold)
        ]
        return self.filtered_viruses

    def compute_distances(self, out_dir):
        """
        Compute distances based on the matches DataFrame in a symmetrical manner.
        
        :return: DataFrame with computed distances.
        """

        self.matches_df['qseqid_genome_ID'] = self.matches_df['qseqid_genome_ID'].astype(str)
        self.matches_df['sseqid_genome_ID'] = self.matches_df['sseqid_genome_ID'].astype(str)
        self.matches_df['original_qseqid_genome_ID'] = self.matches_df['qseqid_genome_ID']
        self.matches_df['original_sseqid_genome_ID'] = self.matches_df['sseqid_genome_ID']
        
        self.matches_df['genome_pair'] = self.matches_df.apply(
            lambda row: tuple(sorted([row['qseqid_genome_ID'], row['sseqid_genome_ID']])),
            axis=1
        )
        
        # Group by the symmetrical pairs and aggregate data
        df = self.matches_df.groupby('genome_pair').agg({
            'conserved_AA_#': 'sum', 
            'qseqid_genome_length': 'first',  
            'sseqid_genome_length': 'first',
            'qseqid_input_seq': 'first',     
            'sseqid_input_seq': 'first',      
            'sseqid_virus': 'first',
            'length': 'sum',
            'pident': 'sum',
            'original_qseqid_genome_ID': 'first',  
            'original_sseqid_genome_ID': 'first'  
        }).reset_index()

        # Remove self-hits
        df = df[df['genome_pair'].apply(lambda x: x[0] != x[1])]
        
        # Calculate average genome length
        df['avg_genome_length'] = (df['qseqid_genome_length'] + df['sseqid_genome_length']) / 2
        
        # Calculate symmetric distance
        df['distance'] = (1 - (df['conserved_AA_#'] * 3) / df['avg_genome_length']).round(6)
        df['conserved_AA_#'] = df['conserved_AA_#'].abs().round(2)
        
        df[['qseqid_genome_ID', 'sseqid_genome_ID']] = pd.DataFrame(df['genome_pair'].tolist(), index=df.index)
        
        df = df.drop(columns=['genome_pair'])

        df = df.sort_values(by=['qseqid_genome_ID', 'distance'])

        self.distances_df = df
        output_file = os.path.join(out_dir, 'distances.csv')
        self.distances_df.to_csv(output_file, index=False)
        return self.distances_df

    def get_top_hits(self, top_n=5):
        """
        Find the top N hits for each genome based on minimum distance and filter them using percentage_viral criteria.
        
        :param top_n: Number of top hits to retain for each genome.
        :return: Filtered DataFrame with top hits.
        """

        if not hasattr(self, 'distances_df'):
            raise ValueError("Distances have not been computed. Please run compute_distances() first.")
        
        self.filter_viruses()

        top_hits = self.distances_df.groupby('qseqid_genome_ID').head(top_n)
        
        self.top_hits_filtered = pd.merge(top_hits, self.filtered_viruses, 
                                          left_on='sseqid_genome_ID', 
                                          right_on='genome_ID')
        
        self.top_hits_filtered.to_csv(os.path.join(self.out_dir, 'top_hits_viruses.csv'), index=False)
        return self.top_hits_filtered

    def create_phylip_table(self):
        """
        Create a PHYLIP table of distances between selected viruses.
        
        :return: String representation of the PHYLIP table.
        """

        if not hasattr(self, 'top_hits_filtered'):
            raise ValueError("Top hits have not been filtered. Please run get_top_hits() first.")
        
        total_steps = 4
        pbar = tqdm(total=total_steps, desc='Converting distance matrix to PHYLIP')

        matrix = self.distances_df.pivot_table(
            columns='qseqid_genome_ID', 
            index='sseqid_genome_ID', 
            values='distance'
        ).reset_index()
        matrix = matrix.set_index('sseqid_genome_ID')
        pbar.update(1)
        
        # Check for and handle missing rows and columns
        cols = set(matrix.columns.to_list())
        indx = set(matrix.index.to_list())
        cols_missing = list(set(indx) - set(cols))
        idx_missing = list(set(cols) - set(indx))
        
        # Add missing rows
        new_row = pd.Series(np.full(shape=len(idx_missing), fill_value=np.nan), index=idx_missing)
        indexes_insert = pd.DataFrame(index=new_row.index)
        matrix = pd.concat([matrix, indexes_insert], axis=1)
        pbar.update(1)
        # Add missing columns
        missing_df = pd.DataFrame(1, columns=cols_missing, index=matrix.index)
        matrix = pd.concat([matrix, missing_df], axis=1)
        sorted_index = sorted(matrix.index)
        matrix = matrix.loc[sorted_index]
        matrix = matrix[sorted_index]
        matrix.columns = [''] * len(matrix.columns)
        matrix.index.names = [len(matrix)]
        matrix.index = matrix.index.str.replace(' ', '_')
        matrix.index = matrix.index.str.replace(r'\,|\(|\)|\:', '_', regex=True)
        pbar.update(1)
        
        matrix.index = matrix.index.str[:64]
        np.fill_diagonal(matrix.values, 0)
        matrix = matrix.fillna(1).astype(float).round(5).clip(lower=0)
        
        phylip_file = os.path.join(self.out_dir, 'phylip.dist')
        matrix.to_csv(phylip_file, header=True, index=True, sep=' ')
        
        pbar.update(1)
        pbar.close()
        print('PHYLIP file saved successfully.')
        return phylip_file


def parse_arguments():
    parser = argparse.ArgumentParser(description='Processes metadata and matches, then calculates viral percentages.')
    parser.add_argument('-p', '--path_to_prots', required=True, help='Path to protein fasta files')
    parser.add_argument('-o', '--output_dir', required=True, help='Directory to save output files')
    parser.add_argument('-m', '--matches', required=True, help='Path to matches.tsv from Diamond')
    parser.add_argument('-md', '--metadata', required=True, help='Path to database metadata')
    return parser.parse_args()


def main():
    MARKDOWN = """
                                 ____                    _  __             _____
                                / ___| _ __ __ _ ___ ___(_)/ _|_   _      /̲\̲ ̲_̲ ̲/̲\ 
                                | |   | '__/ _` / __/ __| | |_| | | |    |__/ \__|
                                | |___| | | (_| \__ \__ \ |  _| |_| |     \/___\/
                                \____ |_|  \__,_|___/___/_|_|  \__, |       |_|
                                                                |___/      /   \\
    Viral taxonomy tool that finds relatedness between viral metagenomes by summing total pairwise protein alignments.
    """
    md = Markdown(MARKDOWN)
    console = Console()
    console.print(md)

    parser = argparse.ArgumentParser()
    args = parse_arguments()

    genome_proteome_mapper = GenomeProteomeMapper(args.path_to_prots, args.output_dir)
    generated_metadata = genome_proteome_mapper.map_genome_to_proteome()
    merger = MetadataMerger(matches=args.matches, generated_metadata=generated_metadata, output_dir=args.output_dir)
    merger.load_and_concat_metadata(input_metadata_path=args.metadata)
    merged_matches = merger.merge_metadata_matches()
    viral_percentage_calculator = ViralPercentageCalculator(merged_matches, generated_metadata, args.output_dir)
    percentage_viral_df = viral_percentage_calculator.calculate_percentage_viral()

    calculator = VirusDistanceCalculator(merged_matches, percentage_viral_df, args.output_dir)
    calculator.compute_distances(args.output_dir)
    calculator.get_top_hits(top_n=5)
    calculator.create_phylip_table()

if __name__ == '__main__':
    main()
