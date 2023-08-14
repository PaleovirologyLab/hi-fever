import argparse

def stop_handler(sequence, action):
    codon_table = {
        "TAA": "taa",
        "TAG": "tag",
        "TGA": "tga"
    }

    stop_codons = 0
    result = []
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i + 3]
        if codon in codon_table:
            stop_codons += 1
            if action == "convert":
                result.append(codon_table[codon])
        else:
            result.append(codon)

    return "".join(result), stop_codons

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Finds in-frame stops in Genewise coding sequence output, and either converts to lower-case or removes them, counting instances.")
    parser.add_argument("--task", choices=["convert", "remove"], required=True, help="Task to perform on in-frame stop codons: 'convert' to convert them to lower case, 'remove' to delete them.")
    parser.add_argument("--file", required=True, help="Path to the input text file.")

    args = parser.parse_args()

    with open(args.file, "r") as file:
        for line in file:
            columns = line.strip().split("\t")
            coding_sequence = columns[13]  # Column 14 (0-based index)
            modified_sequence, stop_codons = stop_handler(coding_sequence, args.task)
            columns[13] = modified_sequence  # Replace the original sequence with the modified one
            print("\t".join(columns) + "\t" + str(stop_codons))
