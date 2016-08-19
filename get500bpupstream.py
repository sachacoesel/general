# Ryan Groussman
# 11 May 2016
# Purpose: Given a list of gene IDs for an organism (initially hard-coded for T. pseudonana),
# use the GFF3 file for these gene IDs to reference the chromosome contigs and retrieve
# the 500bp region upstream of the start codon. Outputs a FASTA file with upstream region for each given gene ID.

# Will capture introns unlike retrieving gene models.

# Usage: get500bpupstream.py inlist outfile.fasta


## Hard-coded filenames for script testing:
# List of gene IDs
list_name = "/Users/rgroussman/Dropbox/Armbrust/bioinfo/scripts/get500upstream_test/ccm-pr_list.txt"

# import packages
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys

only_upstream = False
list_name = sys.argv[1]
outfilename = sys.argv[2]
if len(sys.argv) > 3:
    only_upstream = sys.argv[3]
    # If only_500upstream = True, will only output 500bp upstream without the remaining gene model


## Hard-coded values for GFF, rubric, and chromosome contigs ##
# chromosome name matching rubric, matches GFF and FASTA chromosome names
"/Users/rgroussman/Dropbox/Armbrust/bioinfo/scripts/get500upstream_test/chr_rubric.tab" # not used; hard-coded in script.
# GFF file
gff_name = "/Users/rgroussman/Dropbox/Armbrust/share/data/seq/organisms/Thalassiosira_pseudonana/genome/Thaps3_geneModels_FilteredModels2.gff3"
# chromosome contigs
chr_contigs = "/Users/rgroussman/Dropbox/Armbrust/share/data/seq/organisms/Thalassiosira_pseudonana/genome/Thaps3_assembly_chromosomes.fasta"


def isolate_gene_id(attributes):
    """Given a GFF attributes string specific to Thaps, returns just the
    numerical gene id. e.g, given "ID=867;Name=fgenesh1_pg.C_chr_1000001", returns "867"."""

    att_list = attributes.split(";")
    gene_id = att_list[0]
    gene_id = gene_id[3:]
    return gene_id

def matrix_of_gene_coordinates(gff_name):
    """Given a GFF filename, access the gff and and creates a dictionary with
    the chromosome, start, end, and strand for each gene model in the GFF file.
    Returns a dictionary with these values."""

    # Go through each line of the GFF file and split it by tab. If the feature field is mRNA, then record the start and end positions in a dictionary and strand.

    gene_coordinates_dict = {}
    gff_file = open(gff_name, "r")
    gff_file.readline() # Skip the header
    for line in gff_file:
        line_elts = line.split("\t")
        if line_elts[2] == "mRNA":
            gene_id = isolate_gene_id(line_elts[8])
            gene_coordinates_dict[gene_id] = {}
            gene_coordinates_dict[gene_id]["Chromosome"] = line_elts[0]
            gene_coordinates_dict[gene_id]["Start"] = line_elts[3]
            gene_coordinates_dict[gene_id]["End"] = line_elts[4]
            gene_coordinates_dict[gene_id]["Strand"] = line_elts[6]

    return gene_coordinates_dict

def get_coordinates_from_gene_list(gene_coordinates_dict, gene_list):
    """Given a gene_coordinates_dict and a list of gene IDs, return a smaller
    dictionary with the start, end, and strand data for genes on the list"""

    gene_list_dict = {}
    for gene in gene_list:
        gene = gene.strip()
        if gene in gene_coordinates_dict.keys():
            gene_list_dict[gene] = gene_coordinates_dict[gene]
    return gene_list_dict

def print_gene_list_coords(gene_list_dict):
    """From a given gene_list_dict, will print out a summary table of the
    start, end, strand in a readable format."""

    for gene in gene_list_dict:
            print "###", gene, "###"
            print "Start:", gene_list_dict[gene]["Start"]
            print "End:", gene_list_dict[gene]["End"]
            print "Chromosome:", gene_list_dict[gene]["Chromosome"]
            print "Strand:", gene_list_dict[gene]["Strand"]
            print

def build_fasta_to_gff_chr_names():
    """Builds and returns the dictionaries that link the FASTA and GFF
    variations for chromosome names. Needs only be called once per run."""

    # First element of tuple is FASTA format, second element is GFF format
    chr_rubric = [("chr_1", "Chr1"),
    ("chr_2", "Chr2"),
    ("chr_3", "Chr3"),
    ("chr_4", "Chr4"),
    ("chr_5", "Chr5"),
    ("chr_6", "Chr6"),
    ("chr_7", "Chr7"),
    ("chr_8", "Chr8"),
    ("chr_9", "Chr9"),
    ("chr_10", "Chr10"),
    ("chr_11a", "Chr11a"),
    ("chr_11b", "Chr11b"),
    ("chr_12", "Chr12"),
    ("chr_13", "Chr13"),
    ("chr_14", "Chr14"),
    ("chr_15", "Chr15"),
    ("chr_16a", "Chr16a"),
    ("chr_16b", "Chr16b"),
    ("chr_17", "Chr17"),
    ("chr_18", "Chr18"),
    ("chr_19a_19", "Chr19a_19"),
    ("chr_19b_31", "Chr19b_31"),
    ("chr_19c_29", "Chr19c_29"),
    ("chr_20", "Chr20"),
    ("chr_22", "Chr22"),
    ("chr_23", "Chr23"),
    ("chr_24", "Chr24")]

    # Build dictionaries from the rubric for going either direction
    fasta_to_gff_dict = {}
    for (fasta, gff) in chr_rubric:
        fasta_to_gff_dict[fasta] = [gff]

    gff_to_fasta_dict = {}
    for (fasta, gff) in chr_rubric:
        gff_to_fasta_dict[gff] = [fasta]

    return fasta_to_gff_dict, gff_to_fasta_dict

def fasta_to_gff_chr_names(dict, name):
    """Will convert from FASTA to GFF names for chromosomes; they use slightly different
    syntax. Give the provided dict ('fasta_to_gff_dict' or 'gff_to_fasta_dict') and chr_name, will return the chr_name for the opposite type.
    build_fasta_to_gff_chr_names must have already been called"""

    return dict[name][0]

def first_500(seq):
    """Given a sequence, will return only the first 500 characters. Intended to be used
    when script called with only_upstream=True and only upstream region output desired."""

    return seq[:500]


def capture_flanking_coordinates(gene_list_dict, chr_contigs):
    """Given a dictionary of coordinates from a gene list and a large FASTA
    file with chromosomal contigs, will write out a FASTA file with flanking
    genomic coordinates hard-coded below"""

    # If gene is on - strand, START from higher number and move down to end.
    # Therefore, UPSTREAM is +500 from highest end to start.
    # If gene is on + strand, start from LOWEST number and move up to end.
    # Therefore, UPSTREAM is -500 from lowest end to start.

    # For testing purposes, will print out 500bp upstream + whole CDS in order to
    # compare to traditional sequence retrieval techniques. Can chop off later
    # if only interested in upstream regions.

    # Because the chromosomal sequences are large, we will iterate through each one
    # in turn, and then search for the particular genes in gene_list_dict that need
    # to reference the given contig.


    output_handle = open(outfilename, "w")

    START_BUFFER = 500 # Number of bp to capture before transcription start
    END_BUFFER = 0 # Number of bp to capture after transcription end

    for seq_record in SeqIO.parse(chr_contigs, "fasta"):
        fasta_chr = seq_record.description
        fasta_to_gff_dict, gff_to_fasta_dict = build_fasta_to_gff_chr_names()
        gff_chr = fasta_to_gff_chr_names(fasta_to_gff_dict, seq_record.description)
        for gene in gene_list_dict:
            if gene_list_dict[gene]["Chromosome"] == gff_chr:
                # Create  new SeqRecord
                out_handle = gene
                if gene_list_dict[gene]["Strand"] == "+":
                    start_cut = int(gene_list_dict[gene]["Start"]) - 1 - START_BUFFER
                    end_cut = int(gene_list_dict[gene]["End"]) + END_BUFFER
                    out_sequence = seq_record.seq[start_cut:end_cut]
                elif gene_list_dict[gene]["Strand"] == "-":
                    start_cut = int(gene_list_dict[gene]["Start"]) - 1 - END_BUFFER
                    end_cut = int(gene_list_dict[gene]["End"]) + START_BUFFER
                    out_sequence = seq_record.seq[start_cut:end_cut]
                    out_sequence = out_sequence.reverse_complement()
                if only_upstream == "True":
                    out_sequence = out_sequence[:500]
                out_record = SeqRecord(out_sequence, id= out_handle)
                SeqIO.write(out_record, output_handle, "fasta")
    output_handle.close()

def main():
    # Assemble gene coordinates from GFF
    gene_coordinates_dict = matrix_of_gene_coordinates(gff_name)
    gene_list = open(list_name, "r")

    # Gather coordinates from small list of genes
    gene_list_dict = get_coordinates_from_gene_list(gene_coordinates_dict, gene_list)

    # print_gene_list_coords(gene_list_dict)

    # Build linking dictionaries, capture desired genomic coordinates
    fasta_to_gff_dict, gff_to_fasta_dict = build_fasta_to_gff_chr_names()
    capture_flanking_coordinates(gene_list_dict, chr_contigs)


if __name__ == "__main__":
	main()
