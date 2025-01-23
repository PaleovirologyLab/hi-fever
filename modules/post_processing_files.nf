process merge_seqs_loci {

    input:
    path strict_fastas_collected
    path context_fastas_collected

    output:
    path "loci-merged-coordinates.fasta", emit: loci_merged_fa
    // path "loci-merged-coordinates.fasta.gz", emit: loci_merged_fa_gz
    path "loci-context-coordinates.fasta.gz", emit: loci_merged_context_gz
    path "all_context_coords.bed", emit: all_context_coords_bed

    // publishDir "${params.outdir}", mode: "copy", pattern: "loci-merged-coordinates.fasta"
    // publishDir "${params.outdir}", mode: "copy", pattern: "all_context_coords.bed"

    """

    # Concatenate loci extracted from strictly overlapping coordinates

    cat $strict_fastas_collected > loci-merged-coordinates.fasta

    # Concatenate loci extracted from contextual coordinates
    # Simultaneously, generate accurate context coordinates in BED format 
    # for downstream genewise (since upstream awk operation can leave non-existent overhang on the right flank)

    cat $context_fastas_collected | \
    tee >(grep ">" | sed 's/ .*//; s/>//; s/[:-]/\t/g' | sort -k1,1 -k2,2n > all_context_coords.bed) | \
    gzip > \
    loci-context-coordinates.fasta.gz

    # Zip up query file for outdir publication
    # gzip loci-merged-coordinates.fasta > loci-merged-coordinates.fasta.gz

    """

}

process concatenate_publish_tables {
    input:
    path collected_files
    val table_name  // Receive the table name as a string

    output:
    path "${table_name}"
    publishDir "${params.outdir}/sql", mode: "copy", pattern: "${table_name}"

    """
    cat $collected_files > '${table_name}'
    """
}

// Old way of publishing predicted_orfs and locus_assembly_map
// process publish {

//     input:
//     path locus_assembly_map_collected
//     path orfs_collected

//     output:
//     path "locus_assembly_map.tsv"
//     path "predicted_ORFs.tsv"
//     publishDir "${params.outdir}/sql", mode: "move", pattern: "locus_assembly_map.tsv"
//     publishDir "${params.outdir}/sql", mode: "move", pattern: "predicted_ORFs.tsv"

//     """

//     cat $locus_assembly_map_collected > locus_assembly_map.tsv
//     cat $orfs_collected > predicted_ORFs.tsv

//     """

// }
