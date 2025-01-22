#!/usr/bin/env python

description = """
Retrieve taxonomical information from protein accession numbers
"""

from Bio import SeqIO, Entrez
from glob import glob
from time import sleep

import modin.pandas as pd
import argparse

def get_record_taxonomy(record):
    try:
        tax_info = record.annotations['taxonomy']
        return tax_info
    except Exception as e:
        print(f"Error fetching taxonomy for TaxID {record.id}: {e}")
        return None


if __name__ =="__main__":
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('diamond_hits', type=argparse.FileType('r'),
                        help='table with concatenated hits from forward and reciprocal diamonds: mixed_hits.txt')
    parser.add_argument('email', help='email for Entrez API transactions', type=str)
    parser.add_argument("--batch", type=int, default=50,
                        help="Batch size for retrieving records")
    parser.add_argument("--delay", type=float, default=1.0, 
                        help="Number of seconds to pause between queries (default 1).")
    parser.add_argument("--outfile", type=str, default="diamond-matches-taxonomy.csv",
                        help="Name of the output file")
    args = parser.parse_args()

    # Read diamond hits

    all_hits = pd.read_csv(args.diamond_hits, sep="\t", usecols=[0,1], header=None)
    all_hits.columns = ["hit_accn", "locus"]
    unique_accn = all_hits["hit_accn"].unique()
    accns = unique_accn.tolist()

    print(accns)
    # Get taxonomical information for unique accessions
    Entrez.email = "lmuoz@uwo.ca"
    accns_tax = {}
    for i in range(0, len(accns), args.batch):
        query = ','.join(accns[i:(i+args.batch)])
        response = Entrez.efetch(db="protein", rettype="gb", retmode="text", id=query)
        
        for record in SeqIO.parse(response, format="genbank"):
            taxonomy = get_record_taxonomy(record)
            
            if taxonomy is None: 
                continue
            accns_tax.update({record.id:taxonomy})
        sleep(args.delay)

    final_taxonomy = pd.DataFrame(accns_tax)

    df = pd.DataFrame.from_dict(accns_tax, orient='index')

    # Reset the index to turn keys into a column
    df.reset_index(inplace=True)

    # Save to a CSV file
    df.to_csv(args.outfile, index=False, sep="\t")
