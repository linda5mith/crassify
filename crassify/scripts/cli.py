import argparse
from rich.console import Console
from rich.markdown import Markdown
from crassify.scripts import GenomeProteomeMapper, MetadataMerger, ViralDistanceCalculator, ViralPercentageCalculator, make_summary_plot
import os

def main():
    parser = argparse.ArgumentParser(description="Crassify: viral taxonomy tool")
    parser.add_argument("-p", "--path_to_prots", required=True)
    parser.add_argument("-o", "--output_dir", required=True)
    parser.add_argument("-m", "--matches", required=True)
    parser.add_argument("-md", "--metadata", required=True)
    args = parser.parse_args()

    MARKDOWN = """
                                 ____                    _  __             _____
                                / ___| _ __ __ _ ___ ___(_)/ _|_   _      /̲\̲ ̲_̲ ̲/̲\ 
                                | |   | '__/ _` / __/ __| | |_| | | |    |__/ \__|
                                | |___| | | (_| \__ \__ \ |  _| |_| |     \/___\/
                                \____ |_|  \__,_|___/___/_|_|  \__, |       |_|
                                                                |___/      /   \\
    A viral classification tool that assigns taxonomy to metagenomic contigs by comparing entire proteomes to exemplar viral genomes, measuring relatedness through summed pairwise protein alignments.
    """
    md = Markdown(MARKDOWN)
    console = Console()
    console.print(md)

    gpm = GenomeProteomeMapper(args.path_to_prots, args.output_dir)
    generated_metadata = gpm.map_genome_to_proteome()

    merger = MetadataMerger(args.matches, generated_metadata, args.output_dir)
    merger.load_and_concat_metadata(args.metadata)
    merged_matches = merger.merge_metadata_matches()
    
    vpc = ViralPercentageCalculator(merged_matches, generated_metadata, args.output_dir)
    percentage_viral_df = vpc.calculate_percentage_viral()

    vdc = ViralDistanceCalculator(merged_matches, percentage_viral_df, args.output_dir)
    vdc.compute_distances()
    #vdc.create_phylip_table()


if __name__ == "__main__":
    main()
