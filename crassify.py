from Bio import SeqIO
import pandas as pd
import numpy as np
import os
import re
import argparse
import subprocess
from io import StringIO
from rich.console import Console
from rich.progress import track, Progress
from datetime import datetime

def map_genome_to_proteome(path_to_prots):
    '''Takes a folder of protein fasta files as input, 
    maps the genome to protein ID's and gets the genome length.'''
    if os.path.isdir(path_to_prots):
        prots = list_files(path_to_prots)
        if prots:
            console = Console()
            console.log('Mapping genome ID to protein ID\'s.', style="bold green")
            proteins={}
            genome_lens={}
            for seq in track(prots, description='Parsing proteins...'):
                nt_ID = get_basename(seq)
                records = list(SeqIO.parse(seq,'fasta'))
                record_lens = [len(record.seq) for record in records]
                genome_len = sum(record_lens)*3
                record_ids = [record.id for record in records]
                proteins[nt_ID] = record_ids
                genome_lens[nt_ID] = genome_len
            genome_protein_map = pd.DataFrame.from_dict(proteins, orient='index')
            genome_lens = pd.DataFrame.from_dict(genome_lens, orient='index')
            genome_lens.columns = ['genome_length']
            protein_metadata = pd.merge(genome_protein_map, genome_lens, left_index=True, right_index=True)
            protein_metadata = protein_metadata.reset_index()
            protein_metadata = protein_metadata.rename(columns={'index': 'ID'})
            # Transpose df to have columns ncbi_acc, prot_acc, genome_length
            melt = pd.melt(protein_metadata, id_vars=['ID','genome_length'])
            melt = melt[melt.value.notna()]
            melt = melt[['ID','genome_length','value']] \
                .sort_values(by='ID') \
                .reset_index(drop = True)
            return melt
    else:
        console.print(f'{path_to_prots} is empty or does not contain fasta files.\nCheck -p path/to/proteins.', style="bold underline red")

def merge_metadata_matches(matches, input_metadata, db_metadata):
    '''Merges protein metadata with diamond output matches.tsv'''
    console = Console()
    console.log('Merging metadata with matches.tsv protein alignments.',style="bold green")
    with Progress() as pb:
        t1 = pb.add_task(description='Merging...', total=100)
        matches = pd.read_csv(matches, sep='\t', header=None, low_memory=False)
        matches.columns = ['qseqid','sseqid', "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
        matches = matches.merge(input_metadata, left_on=['qseqid'], right_on=['value']) \
            .rename(columns={'genome_length':'qseqid_genome_length','ID':'qseqid_ID'}) \
            .drop('value', axis=1)
        pb.update(task_id=t1, advance=50)
        db_metadata = pd.read_csv(db_metadata, low_memory=False)
        matches = matches.merge(db_metadata, left_on=['sseqid'], right_on=['protein_ncbi_acc'], how='left') \
            .rename(columns={'genome_length':'sseqid_genome_length'}) 
        pb.update(task_id=t1, advance=50)
    return matches

def calc_dist(df):
    '''Finds best hits (most conserved AA's shared) for each protein in proteome.
    Finds top 5 highest taxonomic matches based on summation of conserved AA's.
    Returns distance matrix based on formula 1-(ORFs_A)+(ORFs_B)/(len(Genome_A+Genome_B)/2)
    '''
    console = Console()
    console.log('Calculating distances between genomes.',style="bold green")
    now = datetime.now()
    dt_string = now.strftime("%d.%m.%Y_%H:%M")
    outpath = f'./crassify_{dt_string}'
    try:
        os.mkdir(outpath)
    except Exception as e:
        print(e)
    with Progress() as pb:
        t1 = pb.add_task('Finding hits...', total=100)
        df['conserved_AA_#'] = (df['length']/100)*df['pident']
        # taxonomy for each protein hit
        df_proteins = df.sort_values(by=['qseqid','sseqid','sstart','send']) \
            .loc[df.reset_index().groupby(['qseqid'])['conserved_AA_#'].idxmax()] \
            .sort_values(by=['qseqid','conserved_AA_#']) \
            .reset_index(drop=True) \
            .to_csv(f'{outpath}/protein_hits.csv',index=False)
        pb.update(task_id=t1, advance=33.33)
        # return top 5 genome hits 
        df_genome_hits = df.groupby(['qseqid_ID','genome_ncbi_acc','virus'])[['conserved_AA_#','length']].sum().reset_index() \
            .sort_values(by=['conserved_AA_#'], ascending=False) \
            .groupby('qseqid_ID').head(5).reset_index(drop=True) \
            .merge(df, left_on=['genome_ncbi_acc'], right_on=['genome_ncbi_acc'], how='left', suffixes=('', '_y')) \
            .drop_duplicates(keep='first').reset_index(drop=True) \
            .groupby(['qseqid_ID','genome_ncbi_acc'])[['definition','qseqid','sseqid','conserved_AA_#_y','pident','length_y']].agg(lambda x: list(x)).reset_index() \
            .merge(df, left_on=['genome_ncbi_acc'], right_on=['genome_ncbi_acc'], how='left', suffixes=('', '_y')) \
            .drop_duplicates(subset=['qseqid_ID','genome_ncbi_acc'], keep='first').reset_index(drop=True)
        df_genome_hits = df_genome_hits.drop(df_genome_hits.filter(regex='_y$').columns, axis=1) \
            .drop('conserved_AA_#', axis=1) \
            .to_csv(f'{outpath}/genome_hits.csv',index=False)
        pb.update(task_id=t1, advance=33.33)
        # calculate distances 
        df_dist = df.groupby(['qseqid_ID','genome_ncbi_acc']).agg({'conserved_AA_#':'sum', 
                                    'qseqid_genome_length':'first',
                                    'sseqid_genome_length':'first',
                                    'virus':'first',
                                    'length':'sum'}).reset_index()
        df_dist['distance'] = 1-(df_dist['conserved_AA_#']*3)/((df_dist['qseqid_genome_length']+df_dist['sseqid_genome_length'])/2)
        df_dist = df_dist.sort_values(by=['qseqid_ID','distance']) \
            .reset_index(drop=True) \
            .round(6) \
            .to_csv(f'{outpath}/dists.csv',index=False)
        pb.update(task_id=t1, advance=34)
    return df_dist 

def to_phylip(df_dist):
    '''Converts distances to PHYLIP format.'''
    console = Console()
    console.print(f'{datetime.now()}: Converting distance matrix to PHYLIP format.')
    matrix = df_dist.pivot_table(columns='qseqid_ID', index='virus', values='distance').reset_index()
    matrix = matrix.set_index('virus')
    #matrix.to_csv('dist_matrix_phylim_raw.dist', header=True, index=True, sep =' ')
    # Check dimensions of matrix
    cols = dist_piv.columns.to_list()
    indx = dist_piv.index.to_list()
    cols_missing = list(set(indx) - set(cols))
    idx_missing = list(set(cols) - set(indx))
    pbar.update(50)
    if len(cols)>len(indx):
        # insert extra column values into index
        for col in cols_missing:
            matrix = matrix.reset_index()
            row_to_append = (matrix.shape[1]-1) * [1]
            row_to_append.insert(0, col)
            matrix.loc[matrix.shape[0]] = row_to_append
            matrix = matrix.set_index('sseqid_genome')
    if len(indx)>len(cols):
        # insert extra index values as columns
        for idx in idx_missing:
            matrix[idx] = 1
    # matrix = matrix.set_index('sseqid_genome')
    #sort dataframe by index and column values
    print(matrix.head())
    sorted_index = sorted(matrix.index)
    matrix = matrix.loc[sorted_index]
    matrix = matrix[sorted_index]
    matrix.columns = [''] * len(matrix.columns)
    matrix.index.names = [len(matrix)]
    # (),: - not permitted by fastme in index
    matrix.index = matrix.index.str.replace('\,|\(|\)|\:','_', regex=True)
    # fastme: taxa names must not be longer than 64 characters
    matrix.index = matrix.index.str[:64]
    np.fill_diagonal(matrix.values, 0)
    matrix = matrix.fillna(1) \
        .astype(float) \
        .round(5) \
        .clip(lower = 0)
    matrix.to_csv('./crassify_output/dist_matrix_phylip.dist', header=True, index=True, sep =' ')

    return matrix

def to_newick(phylip_matrix, path_to_decenttree):
    print(f'{datetime.now()}: Calling DecentTree to convert PHYLIP to Newick format.')
    # ./decenttree -in /data/san/data1/users/linda/crAss_DB/ictv_proposal_seq_translations/dist_matrix_phylim.dist -t NJ -nt 10 -out dist_matrix_NJ.newick
    p = subprocess.call([path_to_decenttree, "-i", phylip_matrix])
    # check p output

def plot_tree(path_to_newick, path_to_metadata):
    # empress tree-plot -t tree.nwk -fm feature-metadata.tsv -o tree-viz
    # feature metadata needs to be tsv for empress
    metadata = read_metadata(path_to_metadata)
    metadata.to_csv('feature-metadata.tsv', sep='\t', index=False)
    return metadata

def is_valid_arg(parser, arg):
    if not os.path.exists(arg):
        parser.error(f'The file {arg} does not exist!')
    else:
        return arg

def list_files(path):
    path_to_files = []
    for root, dir, files in os.walk(path):
        for f in files: 
            if 'fa' in f or 'fasta':
                path_to_files.append(os.path.join(root, f))
    return path_to_files

def get_basename(path_to_file):
    filepath = os.path.splitext(path_to_file)
    base = filepath[0].split('/')[-1]
    return base


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", dest="diamond_matches", required=True,
                        help="path to diamond matches.tsv", metavar="MATCHES")
    parser.add_argument("-p", dest="proteins", required=True,
                        help="path to unclassifed proteomes.", metavar="PROTEINS")
    parser.add_argument("-db", dest="db_metadata", required=True,
                        help="path to database metadata.", metavar="METADATA")
    args = parser.parse_args()

    genome_prot_map = map_genome_to_proteome(args.proteins)
    matches = merge_metadata_matches(args.diamond_matches, genome_prot_map, args.db_metadata)
    dists = calc_dist(matches)
    #phylip = to_phylip(dists)
    # tree = plot_tree('/home/administrator/phd/crass_DB/May_2023_317_proteomes/crassify/dist_matrix_phylim.newick', args.path_to_metadata)



