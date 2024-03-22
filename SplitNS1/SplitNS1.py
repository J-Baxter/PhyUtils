# Divides a fasta alignment into two sub alignments, based on the root of the phylogenetic tree
# NB this is relative, so 'allele A' from one alignment might not necessarily match 'allele A'
# from another.
#
# Arguments: 1) NS alignment
#
# NB: requires Biopython, re & string
#
# Copyright (c) 2024 James Baxter under GNU GENERAL PUBLIC LICENSE Version 3˚=

from Bio import AlignIO
from re import sub
from string import ascii_uppercase
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor


# Read sequence alignment
def parse_args():
    parser = argparse.ArgumentParser(description="a script to do stuff")
    parser.add_argument("alignment")
    args = parser.parse_args()
    return args


# Infer NJ phylogenetic tree using BioPython
def infer_nj(alignment):
    calculator = DistanceCalculator('identity')
    constructor = DistanceTreeConstructor()

    dm = calculator.get_distance(alignment)
    tree = constructor.nj(dm)
    return tree


#  Identify deepest node and extract
def find_clades(tree):
    deepest_node = max(tree.depths().items(), key=lambda x: x[1])[0]

    # Split the alignment
    left_alignment = []
    right_alignment = []

    for clade in tree.find_clades():
        if deepest_node in clade.get_terminals():
            left_alignment.extend(clade.get_terminals())
        else:
            right_alignment.extend(clade.get_terminals())

    return left_alignment, right_alignment


# Write left and right alignments to files
def write_alignments(original_alignment, split_alignments, filepath):
    count = 0
    for split_alignment in split_alignments:
        suffix = '_' + ascii_uppercase[count] + '.fasta'
        file_path = sub('.fasta', suffix, filepath)
        AlignIO.write((original_alignment[seq] for seq in split_alignment), file_path, "fasta")
        count += 1


def main():
    inputs = parse_args()
    input_alignment = AlignIO.read(inputs.alignment, "fasta")
    tree = infer_nj(input_alignment)
    clades = find_clades(tree)
    write_alignments(input_alignment, clades, inputs.alignment)


if __name__ == "__main__":
    main()