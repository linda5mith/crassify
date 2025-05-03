import os
from Bio import Entrez, SeqIO

# Set your email for NCBI Entrez
Entrez.email = "linda.smith@ucc.ie"  # 👈 Replace with your email

# Output folder
output_dir = "test_phages"
os.makedirs(output_dir, exist_ok=True)

# List of famous phages + crAss001
phages = {
    "T4": "NC_000866.4",
    "T7": "NC_001604.1",
    "T3": "NC_003298.1",
    "Lambda": "NC_001416.1",
    "P1": "NC_005856.1",
    "P2": "NC_001895.1",
    "PhiX174": "NC_001422.1",
    "Mu": "NC_000929.1",
    "SP6": "NC_004831.1",
    "N4": "NC_008720.1",
    "MS2": "NC_001417.2",
    "M13": "NC_003287.2",
    "G4": "NC_001420.2",
    "P22": "NC_002371.1",
    "HK97": "NC_002167.1",
    "T5": "NC_005859.1",
    "T1": "NC_005833.1",
    "T2": "NC_003287.1",
    "R17": "NC_001419.1",
    "phiKZ": "NC_004629.1",
    "crAss001": "NC_024711.1"
}

print(f"Downloading {len(phages)} phages to: {output_dir}")

for name, accession in phages.items():
    try:
        print(f"Fetching {name} ({accession})...")
        handle = Entrez.efetch(db="nuccore", id=accession, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        filepath = os.path.join(output_dir, f"{name.replace(' ', '_')}.fna")
        SeqIO.write(record, filepath, "fasta")
        print(f"✅ Saved {name} to {filepath}")
    except Exception as e:
        print(f"❌ Failed to fetch {name}: {e}")
