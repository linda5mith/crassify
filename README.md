# Crassify: Protein-Based Viral Taxonomy Tool

![Crassify Workflow](images/crassify_logo.png)

**Crassify** is a high-throughput, fast tool for computing relatedness between viral genomes using whole-proteome pairwise protein alignments. Designed for metagenomic datasets, it enables:

- Rapid viral species detection
- Novelty detection
- MAG completeness estimation
- Phylogenetic distance estimation

![Crassify Workflow](images/crassify_workflow.png)

Crassify uses **DIAMOND** to align input proteomes against a curated reference database of 14,329 ICTV-classified viral genomes.

---

## Features

- Fast, DIAMOND-based alignment of input proteins to reference DB
- Computes:
  - `% contig viral` — proportion of viral-matching content
  - `% contig completeness` — relative to closest known virus
  - `% proteins aligned` — how much of your proteome was matched
  - Novelty score — penalizes low alignment and completeness
- Detects circularity based on protein start/end overlaps
- Outputs:
  - Pairwise distance matrices
  - PHYLIP-formatted distance tables
  - Tree visualization support via [Empress](https://github.com/biocore/empress)

---

## Input Files

Crassify takes either:
- **Nucleotide sequences** (that get translated), or
- **Protein FASTA files** (`.faa`)  
  (All proteins for a given genome should be supplied together)

---


## 🛠 Installation

### Using Conda/Mamba and `environment.yml`

```bash
git clone https://github.com/linda5mith/crassify.git .
cd crassify
mamba env create -f environment.yml
mamba activate crassify
```

## Running Crassify with Snakemake

After downloading the Crassify reference database, update the paths in your config.yml:
config.yml Example

# Path to DIAMOND reference database
```python
DB: "/absolute/path/to/CRASSIFY_DB.dmnd"

# Metadata file containing protein annotations
metadata: "/absolute/path/to/ICTV_metadata.csv"

# Input files (can be nucleotide or protein)
input_files:
  - path: "sample_data/test_phages_nucl/pooled_test_phages.fna"
    type: "nucleotide"  # or "protein"

# Output directory
output_directory: "/path/to/crassify_output"

# Path to DecentTree binary (for phylogenetic tree generation)
decenttree: "/path/to/decenttree/build/decenttree"
```

# Installing decenttree
```bash
git clone https://github.com/iqbal-lab-org/DecentTree.git
cd DecentTree
mkdir build && cd build
cmake ..
make -j4
```
(Optional) Add to PATH:
```bash
export PATH="$PWD:$PATH"
```
Or specify this path in your config.yml:
```python
decenttree: "/path/to/DecentTree/build/decenttree"
```

## 🧪 Building and Compiling Your Own Reference Database

You can build a custom Crassify-compatible database using your own set of viral proteomes.

### 🔨 Step 1: Create a DIAMOND Database

```bash
diamond makedb --in your_viral_proteomes.faa -d VIRAL_DB.dmnd
```