#!/usr/bin/env python

description = """
Retrieve Lineage information from a TaxID (host lineage)
"""

from Bio import Entrez

import pandas as pd
import argparse

def get_record_lineage(record):
    sci_name = record['ScientificName']
    lineage = record['Lineage']
    taxid = record['TaxId']
    return sci_name, lineage, taxid

if __name__ =="__main__":
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('assembly_table', type=argparse.FileType('r'),
                        help='table with host TaxIDs: assembly_metadata.tsv')
    parser.add_argument('email', help='email for Entrez API transactions', type=str)
    parser.add_argument("--batch", type=int, default=50,
                        help="Batch size for retrieving records")
    parser.add_argument("--delay", type=float, default=1.0, 
                        help="Number of seconds to pause between queries (default 1).")
    parser.add_argument("--outfile", type=str, default="assembly_lineages.tsv",
                        help="Name of the output file")
    args = parser.parse_args()
    
    Entrez.email = args.email
    # Read host taxonomy table
    assembly_md = pd.read_csv(args.assembly_table, sep="\t", header=None)
    assembly_md.columns = [ 'SpeciesName', 'AssemblyAccession', 'FtpPath_GenBank',
                'RefSeq_category', 'AssemblyDescription', 
                'Taxid', 'Organism', 'SpeciesTaxid', 'AssemblyStatus',
                'ContigN50', 'ScaffoldN50', 'Coverage']
    
    # Get unique ids
    taxids = assembly_md['SpeciesTaxid'].unique()
    taxids_list = taxids.tolist()

    # Get taxonomical information for unique accessions
    assembly_lineages = {"SpeciesName":[], "Lineage":[], "TaxId":[]}
    for i in range(0, len(taxids_list), args.batch):
        my_taxids = [str(number) for number in taxids_list[i:(i+args.batch)]]
        query = ','.join(my_taxids)
        response = Entrez.efetch(db="taxonomy", retmode="xml", id=query)
        
        for record in Entrez.read(response):
            sci_name, lineage, taxid = get_record_lineage(record)
            assembly_lineages['SpeciesName'].append(sci_name)
            assembly_lineages['Lineage'].append(lineage)
            assembly_lineages['TaxId'].append(taxid)

    assembly_tax_df = pd.DataFrame.from_dict(assembly_lineages)

    # Save to a CSV file
    assembly_tax_df.to_csv(args.outfile, index=False, sep="\t")
