configfile: "config.yaml"
import yaml
import os
from Bio import SeqIO
from datetime import datetime

# Extract input file type - nucl or prot
input_type = config['input_files'][0]['type']
input_path = config['input_files'][0]['path']
DB_path = config['DB']
DB_metadata = config['metadata']
timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
output_dir = os.path.join(config['output_directory'], f'crassify_out_{timestamp}')
diamond_db_dir = os.path.join(output_dir, 'diamond_db')

print(f"Using DB: {DB_path}")
print(f"Using metadata: {DB_metadata}")

# Create output directory
os.makedirs(output_dir, exist_ok=True)
os.makedirs(diamond_db_dir, exist_ok=True)
print(f"{output_dir} directory created successfully.")

# Validate database
def validate_target_DB(db_path):
    if db_path.endswith('.dmnd'):
        print("Target DB is a DIAMOND database file (.dmnd)")
        return db_path
    elif db_path.endswith('.faa') or db_path.endswith('.fasta'):
        return db_path
    else:
        raise ValueError("Invalid target DB file type. Must be .dmnd, .faa, or .fasta")

validated_DB = validate_target_DB(DB_path)
print("Validated DB:", validated_DB)

# Dynamically resolve input target
def get_target_path():
    if input_type == "nucleotide":
        return os.path.join(diamond_db_dir, 'query_proteomes.faa')
    else:
        return input_path

def get_target(wildcards):
    return get_target_path()

# Master rule
rule all:
    input:
        get_target_path(),
        os.path.join(diamond_db_dir, 'matches.tsv'),
        os.path.join(output_dir, 'dist_matrix_NJ.nwk'),
        os.path.join(output_dir, 'bioempress/tree-viz/')

# Translate nucleotide input to proteins if input is nucleotide
rule translate_metagenome:
    input:
        input_file = input_path
    output:
        translated_proteomes = os.path.join(diamond_db_dir, 'query_proteomes.faa')
    shell:
        "prodigal -i {input.input_file} -a {output.translated_proteomes}"

# DIAMOND blastp
rule run_diamond_blastp:
    input: 
        query_proteins = get_target,
        target_DB = validated_DB
    output:
        matches = os.path.join(diamond_db_dir, 'matches.tsv')
    shell:
        "diamond blastp -d {input.target_DB} -q {input.query_proteins} "
        "--ultra-sensitive --no-self-hits --max-target-seqs 5 "
        "--evalue 0.001 -o {output.matches}"

# Rule for crassify
rule crassify:
    input: 
        matches = os.path.join(diamond_db_dir, 'matches.tsv'),
        input_proteins = get_target,
        target_DB = DB_path,
        metadata = config['metadata']
    output:
        metadata = os.path.join(output_dir, "feature_metadata.tsv"),
        phylip = os.path.join(output_dir, "phylip.dist")
    params:
        config_file = "config.yaml"
    shell:
        """
        python crassify.py \
        -p {input.input_proteins} \
        -m {input.matches} \
        --metadata {DB_metadata} \
        -o {output_dir}
        """

rule decenttree:
    input: 
        phylip=os.path.join(output_dir, 'phylip.dist')
    output:
        newick=os.path.join(output_dir, 'dist_matrix_NJ.nwk')
    shell:
        f"{config['decenttree']} -in {{input}} -t NJ -nt 10 -out {{output}}"

# rule bioempress
rule dendogram:
    input: 
        nwk = os.path.join(output_dir, 'dist_matrix_NJ.nwk'),
        metadata = os.path.join(output_dir,"feature_metadata.tsv")
    output:
        directory(os.path.join(output_dir,"bioempress/tree-viz/"))
    shell:
        "empress tree-plot -t {input.nwk} -fm {input.metadata} -o {output}"

