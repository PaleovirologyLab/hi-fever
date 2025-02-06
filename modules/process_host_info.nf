process download_extract_host_metadata {

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

process get_metadata {
    
    input:
    path assembly_stats

    output:
    path "assembly_metadata.tsv", emit: assembly_metadata
    publishDir "${params.outdir}/sql", mode: "copy", pattern: "assembly_metadata.tsv"

    """
    get_assemblies_metadata.py $assembly_stats $params.email --outfile assembly_metadata.tsv
    
    """

}

process build_host_taxonomy_table {
    
    input:
    path ftp_file
    path assembly_metadata
    path ncbi_taxonomy_table

    output:
    path "host_taxonomy.tsv"
    publishDir "${params.outdir}/sql", mode: "move", pattern: "host_taxonomy.tsv"

    """

    # Extract assembly accessions and remove version number
    sed 's/^https.*\\///g' $ftp_file | sed 's/\\.[0-9]*_.*//g' | sed '/^\$/d' > assembly-ids.txt

    # Grep assemblies being examined and extract species taxid (col. 8 / 7 + 1)    
    grep -Ff assembly-ids.txt $assembly_metadata | awk -F'\\t' '{print \$8}' | sort | uniq > assembly-taxids.txt

    # Find unique taxids and remove N/As
    cat assembly-taxids.txt | sort | uniq > assembly_taxids-uniq.txt
    awk '{if(\$0!~"N/A") print}' assembly_taxids-uniq.txt > assembly_taxid-uniq-nonas.txt

    # Unpack taxonomy dump
    tar -zxvf $ncbi_taxonomy_table

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

process fetch_host_taxonomy {

    input: 
    path assembly_metadata

    output:
    path "host_tax_information.tsv"
    publishDir "${params.outdir}/sql", mode: "copy", pattern: "host_tax_information.tsv"

    """
    get_lineage_from_assembly_id.py $assembly_metadata $params.email --outfile host_tax_information.tsv
    
    """

}

process assembly_stats {

    input:
    path assembly

    output:
    stdout

    """

    stats.sh in=$assembly format=3 addname= | \
    grep -v n_scaffolds | \
    sed 's/\\/.*\\///g; s/_genomic.fna.gz//'

    """

}
process download_assemblies {

        maxForks 10

        input:
        path ftp_dir

        output:
        path "*genomic.fna.gz"

        """

        # Downloads and checks assembly file for corruption, re-attempts if md5 check fails

        max_attempts=5

        count=0

        while read line

                do

                md5check_function () {
                        assemblyFile="\$(echo "\$line" | sed 's/^.*\\///')_genomic.fna.gz"
                        wget -q "\$line/\$assemblyFile" "\$line/md5checksums.txt"
                        grep \$assemblyFile md5checksums.txt > \$assemblyFile.md5; rm md5checksums.txt
                        status=`md5sum -c \$assemblyFile.md5 2>/dev/null | sed 's/.* //'`
                        if [ "\$status" == FAILED ]
                        then
                                if [ "\$count" == "\$max_attempts" ]
                                then
                                        echo "\$assemblyFile FAILED md5check \$max_attempts times, exiting"; exit 1
                                else
                                        echo "\$assemblyFile FAILED md5check, restart function"; rm \$assemblyFile*; ((count=count+1)); md5check_function
                                fi
                        else
                                echo "\$assemblyFile PASSED md5check"; rm \$assemblyFile.md5
                        fi
                }

                md5check_function

                done < $ftp_dir

        """
}

process blastdb {

	// Directives

	debug false
	tag "${meta.id}"

	input:
	tuple val(meta), path(assembly)

	output:
	tuple val(meta), path("*.gz*nsq")

	"""

    gunzip -c ${assembly} | makeblastdb -in - -out ${assembly} -title ${assembly} -dbtype nucl -parse_seqids

	"""

}

process parse_ftp {

    input:
    path ftp_input

    output:
    path "*.ftp.txt"

    """

    while read line
        do
            assemblyID=`echo \$line | sed 's/^.*\\///'`
            echo \$line > \${assemblyID}.ftp.txt
        done < $ftp_input

    """

}
