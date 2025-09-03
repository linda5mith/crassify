configfile: "config.yaml"
import os
from datetime import datetime

timestamp = datetime.now().strftime("%Y-%m-%d")
output_dir = os.path.join(config['outdir'], f'crassify_out_{timestamp}')
diamond_db_dir = os.path.join(output_dir, "diamond_db")
os.makedirs(output_dir, exist_ok=True)
os.makedirs(diamond_db_dir, exist_ok=True)

rule all:
    input:
        os.path.join(output_dir, "crassify_summary.pdf")

rule translate_metagenome:
    input:
        lambda wc: config["input"]["path"]
    output:
        os.path.join(diamond_db_dir, "query_proteomes.faa")
    shell:
        "prodigal -i {input} -a {output} -p meta -m"

rule run_diamond_blastp:
    input:
        query=lambda wc: os.path.join(diamond_db_dir, "query_proteomes.faa") if config["input"]["type"] == "nucleotide" else config["input"]["path"],
        db=config["db"]["dmnd"]
    output:
        os.path.join(diamond_db_dir, "matches.tsv")
    shell:
        "diamond blastp -d {input.db} -q {input.query} "
        "--ultra-sensitive --no-self-hits --max-target-seqs 1 "
        "--evalue 0.01 -o {output}"

rule crassify:
    input:
        matches=os.path.join(diamond_db_dir, "matches.tsv"),
        query=lambda wc: os.path.join(diamond_db_dir, "query_proteomes.faa") if config["input"]["type"] == "nucleotide" else config["input"]["path"],
        metadata=config["db"]["metadata"]
    output:
        os.path.join(output_dir, "percentage_viral.csv"),
        #os.path.join(output_dir, "feature_metadata.tsv"),
        #os.path.join(output_dir, "phylip.dist")
    params:
        outdir=output_dir
    shell:
        "python scripts/cli.py -p {input.query} -m {input.matches} "
        "--metadata {input.metadata} -o {params.outdir}"

rule visualize:
    input:
        os.path.join(output_dir, "percentage_viral.csv")
    output:
        os.path.join(output_dir, "crassify_summary.pdf")
    run:
        from scripts.visualize import make_summary_plot
        make_summary_plot(input[0], output_dir)
