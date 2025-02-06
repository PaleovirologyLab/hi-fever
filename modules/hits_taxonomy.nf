process BUILD_HITS_TAXONOMY_TABLE {

    input:
    path reciprocal_nr_matches
    path reciprocal_rvdb_matches
    path ncbi_taxonomy_table

    output:
    path "hits_taxonomy.tsv"
    publishDir "${params.outdir}/sql", mode: "copy", pattern: "hits_taxonomy.tsv"

    """
    # Collect taxids from reciprocal diamond searches
    cat $reciprocal_nr_matches $reciprocal_rvdb_matches | awk -F'\\t' '{print \$13}' | \
    sed 's/;/\\n/g' | sort | uniq > reciprocal-taxids.txt

    # Find unique taxids and remove N/As
    cat reciprocal-taxids.txt | sort | uniq > all-taxids-uniq.txt
    awk '{if(\$0!~"N/A") print}' all-taxids-uniq.txt > all-taxids-uniq-nonas.txt

    # Unpack taxonomy dump
    tar -zxvf $ncbi_taxonomy_table

    # Use taxonkit to obtain the lineage and print it in tabular format
    taxonkit lineage all-taxids-uniq-nonas.txt --data-dir . | \
    taxonkit reformat --data-dir . -I 1 -f "{k},{K},{p},{c},{o},{f},{g},{s}" -r "N/A" -a | \
    cut -f 1,3 | sed 's/,/\\t/g' > hits_taxonomy.tsv

    # Add N/A row for sequences without classification
    echo \$'N/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A' >> hits_taxonomy.tsv

    # Clean up larger files
    rm *.dmp *.prt taxdump.tar.gz*

    """

}
 
process FETCH_HITS_TAXONOMY_FROM_ACCNS {

    input: 
    path all_diamond_hits

    output:
    path "hits_taxonomy.tsv"
    publishDir "${params.outdir}/sql", mode: "copy", pattern: "hits_taxonomy.tsv"

    """
    get_taxonomy_from_prot_accn.py $all_diamond_hits $params.email --outfile "hits_taxonomy.tsv"
    
    """

}