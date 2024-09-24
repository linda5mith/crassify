configfile: "config.yaml"
import yaml
import os
from Bio import SeqIO

# Extract input file type - nucl or prot
input_type = config['input_files'][0]['type']
input_path = config['input_files'][0]['path']
DB_path = "database/CRASSIFY_DB.dmnd"
DB_metadata = "database/DB_metadata.csv"
output_dir = os.path.join(config['output_directory'], 'crassify_out')

print(input_path)
print(input_type)

# Create output directory
try:
    os.makedirs(output_dir, exist_ok=True)
    print(f"{output_dir} directory created successfully.")
except Exception as e:
    print(e)

def validate_target_DB(db_path):
    if os.path.basename(db_path) == 'ICTV_DB.dmmd':
        print("Target DB is ICTV_DB.dmmd")
        return db_path
    elif db_path.endswith('.dmnd'):
        print("Target DB is a DIAMOND database file (.dmnd)")
        return db_path
    elif is_fasta_file(db_path):
        return db_path
    else:
        print("Invalid target DB file type. Must be .dmnd, .faa, or .fasta")
        raise ValueError("Invalid target DB file type. Must be .dmnd, .faa, or .fasta")

# Set target based on input type
if input_type == "nucleotide":
    target = os.path.join(output_dir, 'diamond_db/query_proteomes.faa')
else:
    target = input_path

print(f"Input TYPE is {input_type}.")
print(f"TARGET {target}.")

if input_type == "nucleotide":
    rule all:
        input:
            os.path.join(output_dir, 'diamond_db/query_proteomes.faa'),  # Translation target
            os.path.join(output_dir, 'diamond_db/matches.tsv'),
            os.path.join(output_dir, 'dist_matrix_NJ.nwk'),
            os.path.join(output_dir, "bioempress/tree-viz/")
else:
    rule all:
        input:
            input_path,  # Protein input
            os.path.join(output_dir, 'diamond_db/matches.tsv'),
            os.path.join(output_dir, 'dist_matrix_NJ.nwk'),
            os.path.join(output_dir, "bioempress/tree-viz/")

# Translate metagenome if input is nucleotide
rule translate_metagenome:
    input:
        input_file=input_path
    output:
        translated_proteomes=os.path.join(output_dir, 'diamond_db/query_proteomes.faa')
    shell:
        "prodigal -i {input.input_file} -a {output.translated_proteomes}"

# Rule for DIAMOND BLASTP
rule run_diamond_blastp:
    input: 
        query_proteins = target,
        target_DB = validate_target_DB(DB_path)
    output:
        matches = os.path.join(output_dir, 'diamond_db/matches.tsv')
    shell:
        "diamond blastp -d {input.target_DB} -q {input.query_proteins} --ultra-sensitive --no-self-hits --max-target-seqs 5 --evalue 0.001 -o {output.matches}"

# Rule for crassify
rule crassify:
    input: 
        matches = os.path.join(output_dir, 'diamond_db/matches.tsv'),
        input_proteins = target,
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

