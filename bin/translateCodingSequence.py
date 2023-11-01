#!/usr/bin/env python

import argparse

# Dictionary for codon to amino acid mapping
codon_table = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TGT": "C", "TGC": "C",
    "TGG": "W", "CTT": "L", "CTC": "L", "CTA": "L",
    "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P",
    "CCG": "P", "CAT": "H", "CAC": "H", "CAA": "Q", 
    "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R",
    "CGG": "R", "ATT": "I", "ATC": "I", "ATA": "I",
    "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T",
    "ACG": "T", "AAT": "N", "AAC": "N", "AAA": "K",
    "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", 
    "AGG": "R", "GTT": "V", "GTC": "V", "GTA": "V",
    "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A",
    "GCG": "A", "GAT": "D", "GAC": "D", "GAA": "E",
    "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G", 
    "GGG": "G"
}

def translate_dna_to_peptide(dna_sequence):
    peptide_sequence = ""
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i + 3]
        amino_acid = codon_table.get(codon, "X")  # Default to "X" if codon is not found
        peptide_sequence += amino_acid
    return peptide_sequence

def process_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            columns = line.strip().split('\t')
            coding_sequence = columns[10]  # Column 11 (0-based index)
            peptide_sequence = translate_dna_to_peptide(coding_sequence)
            # Insert the peptide sequence as the new column 12
            columns.insert(11, peptide_sequence)
            outfile.write("\t".join(columns) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Translate coding DNA sequences to peptides and insert as a new column.")
    parser.add_argument("--input", required=True, help="Path to the input file.")
    parser.add_argument("--output", required=True, help="Path to the output file.")

    args = parser.parse_args()

    process_file(args.input, args.output)
