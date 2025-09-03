import pandas as pd
import matplotlib.pyplot as plt
import os


def make_summary_plot(percentage_csv, output_dir):
    """
    Generate a 4-panel summary PDF for Crassify results.

    Parameters
    ----------
    percentage_csv : str
        Path to percentage_viral.csv (Crassify output).
    output_dir : str
        Directory to save the PDF summary.
    """

    df = pd.read_csv(percentage_csv)

    # Get top 20 most common species hits
    top_hits = df["top_species_hit"].value_counts().head(20)

    # Protein hits raw hist
    protein_hits = pd.to_numeric(df["protein_hits"], errors="coerce")

    # Protein hits binned
    protein_hits = protein_hits.fillna(0)
    bins = [-0.1, 0, 1, 2, 5, 10, protein_hits.max()]
    labels = ["0", "1", "2", "3–5", "6–10", ">10"]
    protein_bins = pd.cut(protein_hits, bins=bins, labels=labels, right=True, include_lowest=True)
    bin_counts = protein_bins.value_counts().reindex(labels, fill_value=0)

    # Contig completeness
    contig_completeness = pd.to_numeric(df["% contig completeness"], errors="coerce").dropna()

    # --- 4-panel plot ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # 1. Top hits barh
    top_hits.plot(kind="barh", color="teal", width=0.8, ax=axes[0, 0])
    axes[0, 0].invert_yaxis()
    axes[0, 0].set_xlabel("Contig Count", fontsize=12)
    axes[0, 0].set_ylabel("Top Virus Hit", fontsize=12)
    axes[0, 0].set_title("Top 20 Viral Species Hits", fontsize=14, loc="left")

    # 2. Raw protein hits hist
    axes[0, 1].hist(protein_hits, bins=50, color="teal", edgecolor="black")
    axes[0, 1].set_xlabel("Viral Protein Hits per Contig", fontsize=12)
    axes[0, 1].set_ylabel("Contig Counts", fontsize=12)
    axes[0, 1].set_title("Distribution of Protein Hits", fontsize=14, loc="left")
    axes[0, 1].set_yscale("log")

    # 3. Protein hits binned barplot
    bin_counts.plot(kind="bar", color="teal", edgecolor="black", ax=axes[1, 0])
    axes[1, 0].set_xlabel("Protein hits", fontsize=12)
    axes[1, 0].set_ylabel("Number of contigs", fontsize=12)
    axes[1, 0].set_title("Protein Hit Counts per Contig", fontsize=14, loc="left")

    # 4. % contig completeness distribution
    axes[1, 1].hist(contig_completeness, bins=50, color="teal", edgecolor="black")
    axes[1, 1].set_xlabel("% Contig Completeness", fontsize=12)
    axes[1, 1].set_ylabel("Number of Contigs", fontsize=12)
    axes[1, 1].set_title("Distribution of Contig Completeness", fontsize=14, loc="left")

    # Add large title
    fig.suptitle(
        "Crassify: Viral Classification Summary for Contigs",
        fontsize=22,
        fontweight="bold",
        y=1.02,
        x=0.505,
    )

    plt.tight_layout()

    # Save
    out_pdf = os.path.join(output_dir, "crassify_summary.pdf")
    plt.savefig(out_pdf, format="pdf", bbox_inches="tight")

    out_png = os.path.join(output_dir, "crassify_summary.png")
    plt.savefig(out_png, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    return out_pdf
