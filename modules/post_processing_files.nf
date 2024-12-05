process find_best_hits {

    input:
    path forward_matches
    path query_proteins
    path reciprocal_matches
    path reciprocal_seqs
    path reciprocal_hits

    output:
    path "matches.dmnd.annot.tsv"
    path "best_hits.fasta", emit: best_hits_fa_ch
    path "best_pairs.txt", emit: best_pairs_txt
    
    publishDir "${params.outdir}", mode: "copy", pattern: "matches.dmnd.annot.tsv"
    publishDir "${params.outdir}", mode: "copy", pattern: "best_hits.fasta"
    publishDir "${params.outdir}", mode: "copy", pattern: "best_pairs.txt"
    publishDir "${params.outdir}", mode: "copy", pattern: "forward_hits.txt"

    // publishDir "${params.outdir}/sql", mode: "copy", pattern: "reciprocal-nr-matches.dmnd.tsv"
    // publishDir "${params.outdir}/sql", mode: "copy", pattern: "reciprocal-rvdb-matches.dmnd.tsv"
    // publishDir "${params.outdir}/accessory_fastas", mode: "copy", pattern: "*.fasta.gz"
    // publishDir "${params.outdir}/sql", mode: "move", pattern: "matches.dmnd.annot.tsv"

    """
    

    cat $forward_matches > matches.dmnd.annot.tsv

    # Extract best forward hit pairs & their alignment length
    awk 'BEGIN{OFS="\t"}; {print \$1,\$4,\$9-\$8, "forward"}' matches.dmnd.annot.tsv | \
    uniq > \
    forward_hits.txt
    
    # Merge reciprocal and forward hits
    cat forward_hits.txt $reciprocal_hits > mixed_hits.txt

    # Determine best single protein-locus pair from both the forward & reciprocal searches

    sort -k2,2 -k3,3nr mixed_hits.txt | \
    sort -u -k2,2 | \
    tee >(grep reciprocal | cut -f1 | sort | uniq > extract_from_reciprocal.txt) | \
    tee >(grep forward | cut -f1 | sort | uniq > extract_from_forward.txt) > \
    best_pairs_partial.txt

    # Add true genomic start for later conversion

    cut -f2 best_pairs_partial.txt | sed 's/[:-]/\t/g' > coords.txt
    paste best_pairs_partial.txt coords.txt > best_pairs.txt

    # Extract best proteins

    seqtk subseq $reciprocal_seqs extract_from_reciprocal.txt | \
    sed 's/ .*//' > \
    best_hits.fasta

    seqtk subseq $query_proteins extract_from_forward.txt | \
    sed 's/ .*//' >> \
    best_hits.fasta

    """

}

process merge_seqs_loci {

    input:
    path strict_fastas_collected
    path context_fastas_collected

    output:
    path "loci-merged-coordinates.fasta", emit: loci_merged_fa
    // path "loci-merged-coordinates.fasta.gz", emit: loci_merged_fa_gz
    path "loci-context-coordinates.fasta.gz", emit: loci_merged_context_gz
    path "all_context_coords.bed", emit: all_context_coords_bed

    publishDir "${params.outdir}", mode: "copy", pattern: "loci-merged-coordinates.fasta"
    publishDir "${params.outdir}", mode: "copy", pattern: "all_context_coords.bed"

    // publishDir "${params.outdir}/sql", mode: "copy", pattern: "reciprocal-nr-matches.dmnd.tsv"
    // publishDir "${params.outdir}/sql", mode: "copy", pattern: "reciprocal-rvdb-matches.dmnd.tsv"
    // publishDir "${params.outdir}/accessory_fastas", mode: "copy", pattern: "*.fasta.gz"
    // publishDir "${params.outdir}/sql", mode: "move", pattern: "matches.dmnd.annot.tsv"

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
