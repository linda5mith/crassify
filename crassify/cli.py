#!/usr/bin/env python3
import argparse, os, yaml, subprocess
from datetime import datetime
from rich.console import Console
from rich.markdown import Markdown
import importlib.resources as ir

def get_default_data_path(fname):
    return str(ir.files("crassify").joinpath("data", fname))

def main():
    parser = argparse.ArgumentParser(description="Crassify: protein-based viral taxonomy tool")
    parser.add_argument("-i", "--input", required=True, help="Path to input FASTA (nucleotide or protein)")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    parser.add_argument("-t", "--type", choices=["nucleotide", "protein"], default="nucleotide")
    parser.add_argument("--db", default=get_default_data_path("crassify_DB_2025.dmnd"),
                    help="Path to DIAMOND DB")
    parser.add_argument("--metadata", default=get_default_data_path("crassify_metadata.csv"),
                        help="Path to DIAMOND DB metadata")
    args = parser.parse_args()

    console = Console()
    console.print(Markdown("# 🧬 Crassify — Protein-based viral taxonomy tool"))

    # Prepare output dirs
    timestamp = datetime.now().strftime("%Y-%m-%d")
    outdir_abs = os.path.abspath(args.outdir)
    os.makedirs(outdir_abs, exist_ok=True)

    # Build config dict
    config = {
        "input": {"type": args.type, "path": os.path.abspath(args.input)},
        "db": {"dmnd": os.path.abspath(args.db), "metadata": os.path.abspath(args.metadata)},
        "outdir": outdir_abs,
        "steps": {"translate": True, "diamond": True, "crassify": True, "visualize": True}
    }

    # Save resolved config
    resolved_config = os.path.join(outdir_abs, "resolved_config.yaml")
    with open(resolved_config, "w") as f:
        yaml.safe_dump(config, f)

    # Run Snakemake
    snakefile_path = os.path.join(os.path.dirname(__file__), "Snakefile")
    cmd = ["snakemake", "--snakefile", snakefile_path, "--configfile", resolved_config, "--cores", "all"]
    console.print("[bold green]Running Snakemake...[/bold green]")
    subprocess.run(cmd, check=True)
