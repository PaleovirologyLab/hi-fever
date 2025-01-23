process build_taxonomy_table {

    input:
    path ftp_file
    path assembly_metadata
    path reciprocal_nr_matches
    path reciprocal_rvdb_matches

    output:
    path "taxonomy_table.tsv"
    publishDir "${params.outdir}/sql", mode: "move", pattern: "taxonomy_table.tsv"

    """

    # Extract assembly accessions and remove version number
    
    sed 's/^https.*\\///g' $ftp_file | sed 's/\\.[0-9]*_.*//g' | sed '/^\$/d' > assembly-ids.txt

    # Grep assemblies being examined and extract species taxid (col. 8 / 7 + 1)
    
    grep -Ff assembly-ids.txt $assembly_metadata | awk -F'\\t' '{print \$8}' | sort | uniq > assembly-taxids.txt

    # Collect taxids from reciprocal diamond searches
    
    cat $reciprocal_nr_matches $reciprocal_rvdb_matches | awk -F'\\t' '{print \$13}' | sed 's/;/\\n/g' | sort | uniq > reciprocal-taxids.txt

    # Find unique taxids and remove N/As
    
    cat reciprocal-taxids.txt assembly-taxids.txt | sort | uniq > all-taxids-uniq.txt
    awk '{if(\$0!~"N/A") print}' all-taxids-uniq.txt > all-taxids-uniq-nonas.txt

    # Get the taxonomy files from NCBI, check for corruption, re-attempt if md5 check fails

    max_attempts=5

    count=0

    md5check_function () {
        wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz.md5
        status=`md5sum -c taxdump.tar.gz.md5 2>/dev/null | sed 's/.* //'`
        if [ "\$status" == FAILED ]
        then
            if [ "\$count" == "\$max_attempts" ]
            then
                echo "Taxonomy FAILED md5check \$max_attempts times, exiting"; exit 1
            else
                echo "Taxonomy FAILED md5check, restart function"; rm taxdump.tar.gz*; ((count=count+1)); md5check_function
            fi
        else
                echo "Taxonomy PASSED md5check"
        fi
    }

    md5check_function

    # Unpack taxonomy dump

    tar -zxvf taxdump.tar.gz

    # Use taxonkit to obtain the lineage and print it in tabular format
    
    taxonkit lineage all-taxids-uniq-nonas.txt --data-dir . | taxonkit reformat --data-dir . -I 1 -f "{k},{K},{p},{c},{o},{f},{g},{s}" -r "N/A" -a | cut -f 1,3 | sed 's/,/\\t/g' > taxonomy_table.tsv

    # Add N/A row for sequences without classification
    
    echo \$'N/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A' >> taxonomy_table.tsv

    # Clean up larger files

    rm *.dmp *.prt taxdump.tar.gz* "\$(readlink -f $reciprocal_nr_matches)" "\$(readlink -f $reciprocal_rvdb_matches)" "\$(readlink -f $assembly_metadata)"

    """

}
 
process build_diamond_taxonomy_table {

    input: 
    path all_diamond_hits

    output:
    path "diamond_taxonomy_table.tsv"
    publishDir "${params.outdir}/sql", mode: "copy", pattern: "diamond_taxonomy_table.tsv"

    """
    get_taxonomy_from_prot_accn.py $all_diamond_hits $params.email --outfile "diamond_taxonomy_table.tsv"
    
    """

}


process build_host_lineage_table {

    input: 
    path assembly_metadata

    output:
    path "assemblies_lineage_information.tsv"
    publishDir "${params.outdir}/sql", mode: "copy", pattern: "assemblies_lineage_information.tsv"

    """
    get_lineage_from_assembly_id.py $assembly_metadata $params.email --outfile assemblies_lineage_information.tsv
    
    """

}