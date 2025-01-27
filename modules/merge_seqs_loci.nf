process merge_seqs_loci {

    input:
    path strict_fastas_collected
    path context_fastas_collected

    output:
    path "loci-merged-coordinates.fasta", emit: loci_merged_fa
    // path "loci-merged-coordinates.fasta.gz", emit: loci_merged_fa_gz
    path "loci-context-coordinates.fasta.gz", emit: loci_merged_context_gz
    path "all_context_coords.bed", emit: all_context_coords_bed

    publishDir "${params.outdir}/accesory_fastas", mode: "copy", pattern: "loci-merged-coordinates.fasta.gz"
    publishDir "${params.outdir}/accesory_fastas", mode: "copy", pattern: "loci-context-coordinates.fasta.gz"

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
    gzip loci-merged-coordinates.fasta > loci-merged-coordinates.fasta.gz

    """

}

