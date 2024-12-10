#!/usr/bin/env python

description = """
Retrieve information and genome metadata from an assemblyID
Modified from a PoonLab script: https://github.com/PoonLab/surfaces/blob/main/scripts/get_all_accns.py
"""

from Bio import Entrez
from time import sleep
from csv import DictReader, DictWriter, reader, writer

import sys
import argparse
import re

def get_all_ids(species_list, batch, delay):
    """
    From a list of species, retrieve to 400 assembly Ids
    
    """
    assembly_ids = []
    for i in species_list:
        handle = Entrez.esearch(db="assembly", term=i, rettype = "ulist", retmax=400)
        record = Entrez.read(handle)
        # print(record)
        idList = record['IdList']
        assembly_ids.extend(idList)
        sleep(1)
    return assembly_ids

def get_relevant_info(record, keys):
    """
    For an assembly record, get the information specified in keys
    """
    record_info = {}
    # Transfor ftp protocol into http link
    for key in keys:

        if key == 'FtpPath_GenBank':
            if record['AssemblyAccession'].startswith('GCF'):
                subject = record['FtpPath_RefSeq'].replace('ftp', 'http',1)
            else: 
                subject = record[key].replace('ftp', 'http',1)     
            record_info[key] = subject
            continue
        
        elif key == 'ExclFromRefSeq':
            reasons = record[key]
            reason = '-'.join(reasons)
            record_info[key] = reason
            continue

        record_info[key] = record[key]
    
    return record_info

def get_assemblies_per_species(ids, batch, info, delay):
    """
    From a list of Ids, return the assembly information
    :para ids: 
    :para batch: how many records to retrieve from Efetch
    :param info: list, keys of relevant information to retrieve from assembly
    :return sps_assemblies: dict, keyed by species, values are dicts with 
    """

    assembly_info = []
  
    for i in range(0, len(ids), batch):
        query = ','.join(ids[i:i+batch])
        handle = Entrez.esummary(db="assembly", id=query, report="full", retmax=500)
        response = Entrez.read(handle)
        records = response['DocumentSummarySet']['DocumentSummary']
        for record in records:    
            record_info = get_relevant_info(record, info)
            assembly_info.append(record_info)

    return assembly_info
    
if __name__ =="__main__":
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('infile', type=argparse.FileType('r'),
                        help='table with accession numbers: assembly_stats.tsv')
    parser.add_argument('email', help='email for Entrez API transactions', type=str)
    parser.add_argument("--batch", type=int, default=50,
                        help="Batch size for retrieving records")
    parser.add_argument("--delay", type=float, default=1.0, 
                        help="Number of seconds to pause between queries (default 1).")
    parser.add_argument("--outfile", type=str, default="assembly_info.csv",
                        help="Name of the output file")
    args = parser.parse_args()

    reader = reader(args.infile, delimiter='\t')

    assembly_id = []
    for row in reader:
        assembly_name = row[19]
        # To do: Find a better regular expression to extract the genome ID (everything befer)
        assembly = assembly_name[0:15]
        assembly_id.append(assembly)

    batch = args.batch
    Entrez.email = args.email
    
    # # Get assembly Ids associated with each species
    # print(f"Getting Id list for {len(sps)} species in {args.infile.name}\n")
    ids = get_all_ids(assembly_id, batch, args.delay)

    # print(f"Getting info for {len(ids)} ids")
    # # Select the columns the relevant information I want to extract from the assembly
    info = [ 'SpeciesName', 'AssemblyAccession', 'FtpPath_GenBank',
            'RefSeq_category', 'AssemblyDescription', 
            'Taxid', 'Organism', 'SpeciesTaxid', 'AssemblyStatus',
            'ContigN50', 'ScaffoldN50', 'Coverage']
    
    # Get information for all assemblies
    all_assemblies = get_assemblies_per_species(ids, args.batch, info, args.delay)

    with open(args.outfile, 'w', newline='') as output_file:
        dict_writer = DictWriter(output_file, info, delimiter=',')
        dict_writer.writerows(all_assemblies)
