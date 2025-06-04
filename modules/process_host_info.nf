process DOWNLOAD_EXTRACT_HOST_METADATA {

    output:
    path "assembly_metadata.tsv", emit: assembly_metadata_ch
    publishDir "${params.outdir}/sql", mode: "copy", pattern: "assembly_metadata.tsv"

    """

    # Download current eukaryotic assembly summaries

    wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/assembly_summary.txt; mv assembly_summary.txt rs_assembly_summary_invert.txt
    wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt; mv assembly_summary.txt rs_assembly_summary_fungi.txt
    wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/assembly_summary.txt; mv assembly_summary.txt rs_assembly_summary_plant.txt
    wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/assembly_summary.txt; mv assembly_summary.txt rs_assembly_summary_protozoa.txt
    wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/assembly_summary.txt; mv assembly_summary.txt rs_assembly_summary_mammal.txt
    wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/assembly_summary.txt; mv assembly_summary.txt rs_assembly_summary_vert_other.txt
    wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/invertebrate/assembly_summary.txt; mv assembly_summary.txt gb_assembly_summary_invert.txt
    wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/assembly_summary.txt; mv assembly_summary.txt gb_assembly_summary_fungi.txt
    wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/assembly_summary.txt; mv assembly_summary.txt gb_assembly_summary_plant.txt
    wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/protozoa/assembly_summary.txt; mv assembly_summary.txt gb_assembly_summary_protozoa.txt
    wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/assembly_summary.txt; mv assembly_summary.txt gb_assembly_summary_mammal.txt
    wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_other/assembly_summary.txt; mv assembly_summary.txt gb_assembly_summary_vert_other.txt

    # Remove header lines and concatenate

    for i in ??_assembly_summary_*; do grep -v "README_assembly_summary\\|#assembly_accession" \$i >> metadata_partial.txt; done

    # Extract AssemblyID

    cut -f20 metadata_partial.txt | sed 's/.*\\///g' > assemblyID

    # Paste

    paste assemblyID metadata_partial.txt > assembly_metadata.tsv

    """

}



process BUILD_HOST_TAXONOMY_TABLE {

	container 'oras://community.wave.seqera.io/library/taxonkit:0.18.0--f6ca10e249e16f7d'
	conda 'bioconda::taxonkit=0.18.0'

    input:
    path ftp_file
    path assembly_metadata
    path ncbi_taxonomy_table

    output:
    path "host_taxonomy.tsv"
    publishDir "${params.outdir}/sql", mode: "copy", pattern: "host_taxonomy.tsv"

    """

    # Extract assembly accessions and remove version number
    sed 's/^https.*\\///g' ${ftp_file} | sed 's/\\.[0-9]*_.*//g' | sed '/^\$/d' > assembly-ids.txt

    # Grep assemblies being examined and extract species taxid (col. 8 / 7 + 1)    
    grep -Ff assembly-ids.txt ${assembly_metadata} | awk -F'\\t' '{print \$8}' | sort | uniq > assembly-taxids.txt

    # Find unique taxids and remove N/As
    cat assembly-taxids.txt | sort | uniq > assembly_taxids-uniq.txt
    awk '{if(\$0!~"N/A") print}' assembly_taxids-uniq.txt > assembly_taxid-uniq-nonas.txt

    # Unpack taxonomy dump
    tar -zxvf ${ncbi_taxonomy_table}

    # Use taxonkit to obtain the lineage and print it in tabular format
    taxonkit lineage assembly_taxid-uniq-nonas.txt --data-dir . | \
    taxonkit reformat --data-dir . -I 1 -f "{k},{K},{p},{c},{o},{f},{g},{s}" -r "N/A" -a | \
    cut -f 1,3 | sed 's/,/\\t/g' > host_taxonomy.tsv

    # Add N/A row for sequences without classification
    echo \$'N/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A' >> host_taxonomy.tsv

    # Clean up larger files
    rm *.dmp *.prt taxdump.tar.gz*
    """

}

