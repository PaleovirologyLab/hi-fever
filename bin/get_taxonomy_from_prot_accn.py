#!/usr/bin/env python

description = """
Retrieve taxonomical information from protein accession numbers
"""

from Bio import SeqIO, Entrez
from glob import glob
from time import sleep

import modin.pandas as pd
import argparse
import sys

def get_record_taxonomy(record):
    try:
        tax_info = record.annotations['taxonomy']
        return tax_info
    except Exception as e:
        print(f"Error fetching taxonomy for TaxID {record.id}: {e}")
        return None

def parse_taxonomy_information (taxonomy):

    tax_info = {}
    tax_info['second_field']=taxonomy[1]
    tax_info['third_field']=taxonomy[2]
    tax_info['fifth_field']=taxonomy[4]
    tax_info['last_field']=taxonomy[-1]
    tax_info['second_to_last']=taxonomy[-2]
    tax_info['third_to_last']=taxonomy[-2]
    tax_info['all_taxonomy']= ':'.join(taxonomy)

    # Get virus information
    tax_info['family']= 'NA'
    tax_info['viral_order']= 'NA'
    tax_info['viral_kingdom']= 'NA'

    for item in taxonomy:
        if item.endswith('dae'):
            tax_info.update({'family': item})
        elif item.endswith('virales'):
            tax_info.update({'viral_order': item})
        elif item.endswith('virae'):
            # print(item, type(item))
            tax_info.update({'viral_kingdom': item})

    return tax_info    


if __name__ =="__main__":
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('diamond_hits', type=argparse.FileType('r'),
                        help='table with concatenated hits from forward and reciprocal diamonds: mixed_hits.txt')
    parser.add_argument('email', help='email for Entrez API transactions', type=str)
    parser.add_argument("--batch", type=int, default=100,
                        help="Batch size for retrieving records")
    parser.add_argument("--delay", type=float, default=1.0, 
                        help="Number of seconds to pause between queries (default 1).")
    parser.add_argument("--outfile", type=str, default="diamond-matches-taxonomy.csv",
                        help="Name of the output file")
    parser.add_argument("--debug", action='store_true', 
                        help="In debug mode take a table of accns as input")
    args = parser.parse_args()

    # Read diamond hits

    if args.debug: 
        all_hits = pd.read_csv(args.diamond_hits, sep="\t", usecols=[0], header=None)
        all_hits.columns = ["hit_accn"]
    else:
        all_hits = pd.read_csv(args.diamond_hits, sep="\t", usecols=[0,1], header=None)
        all_hits.columns = ["hit_accn", "locus"]
    
    unique_accn = all_hits["hit_accn"].unique()
    accns = unique_accn.tolist()
    # Get taxonomical information for unique accessions
    Entrez.email = args.email
    accns_tax = []
    
    for i in range(0, len(accns), args.batch):
        # print(f"batch: {i}")
        query = ','.join(accns[i:(i+args.batch)])
        response = Entrez.efetch(db="protein", rettype="gb", retmode="text", id=query)
    
        for record in SeqIO.parse(response, format="genbank"):
            taxonomy = get_record_taxonomy(record)
                
            if taxonomy is None:
                # Not taxonomical information provided, omit this record 
                continue
            
            tax_info = parse_taxonomy_information(taxonomy)
            tax_info['record_id'] = record.id  # Add record_id to table
            accns_tax.append(tax_info)
            # print(tax_info)
        
        sleep(args.delay)

    df = pd.DataFrame.from_dict(accns_tax)

    # Reset the index to turn keys into a column
    df.reset_index(inplace=True)

    # Save to a CSV file
    df.to_csv(args.outfile, index=False, sep="\t")
