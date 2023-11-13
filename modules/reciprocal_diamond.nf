process reciprocal_diamond {

    input:
    path strict_fastas_collected
    path context_fastas_collected
    path reciprocal_nr_db
    path reciprocal_rvdb_db
    path annotated_hits_collected
    path clustered_proteins

    output:
    path "reciprocal-nr-matches.dmnd.tsv"
    path "reciprocal-rvdb-matches.dmnd.tsv"
    path "best_hits.fasta", emit: best_hits_fa_ch
    path "loci-merged-coordinates.fasta.gz", emit: merged_fa_ch
    path "loci-context-coordinates.fasta.gz", emit: context_fa_ch
    path "all_context_coords.bed", emit: context_coords_ch
    path "matches.dmnd.annot.tsv"
    path "best_pairs.txt", emit: pairs_ch
    publishDir "${params.outdir}/sql", mode: "move", pattern: "reciprocal-nr-matches.dmnd.tsv"
    publishDir "${params.outdir}/sql", mode: "move", pattern: "reciprocal-rvdb-matches.dmnd.tsv"
    publishDir "${params.outdir}/accessory_fastas", mode: "copy", pattern: "*.fasta.gz"
    publishDir "${params.outdir}/sql", mode: "move", pattern: "matches.dmnd.annot.tsv"

    """

    # Concatenate loci extracted from strictly overlapping coordinates

    cat $strict_fastas_collected > loci-merged-coordinates.fasta

    # Concatenate loci extracted from contextual coordinates
    # Simultaneously, generate accurate context coordinates in BED format for downstream genewise (since upstream awk operation can leave non-existent overhang on the right flank)

    cat $context_fastas_collected | \
    tee >(grep ">" | sed 's/ .*//; s/>//; s/[:-]/\t/g' | sort -k1,1 -k2,2n > all_context_coords.bed) | \
    gzip > \
    loci-context-coordinates.fasta.gz

    # Run reciprocal nr search, extract best reciprocal hit pairs & their alignment length

    diamond blastx \
    --$params.diamond_mode \
    --matrix $params.diamond_matrix \
    --masking seg \
    -d $reciprocal_nr_db \
    -q loci-merged-coordinates.fasta \
    -e 1e-5 \
    -k 20 \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms skingdoms sphylums stitle full_sseq | \
    tee >(sort -k1,1 -k12,12nr | sort -u -k1,1 | tee >(awk 'BEGIN{OFS="\t"}; {print \$2, \$1, \$4 * 3, "reciprocal-nr"}' > mixed_hits.txt) | cut -f2,19 | sort -u -k1,1 | awk '{print ">"\$1"\\n"\$2}' > reciprocal_nr_subset.fasta) | \
    cut -f1-18 > \
    reciprocal-nr-matches.dmnd.tsv

    # Run reciprocal RVDB search, extract best reciprocal hit pairs & their alignment length

    diamond blastx \
    --$params.diamond_mode \
    --matrix $params.diamond_matrix \
    --masking seg \
    -d $reciprocal_rvdb_db \
    -q loci-merged-coordinates.fasta \
    -e 1e-5 \
    -k 20 \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms skingdoms sphylums stitle full_sseq | \
    tee >(sort -k1,1 -k12,12nr | sort -u -k1,1 | tee >(awk 'BEGIN{OFS="\t"}; {print \$2, \$1, \$4 * 3, "reciprocal-rvdb"}' >> mixed_hits.txt) | cut -f2,19 | sort -u -k1,1 | awk '{print ">"\$1"\\n"\$2}' > reciprocal_rvdb_subset.fasta) | \
    cut -f1-18 > \
    reciprocal-rvdb-matches.dmnd.tsv

    # Zip up query file for outdir publication

    gzip loci-merged-coordinates.fasta

    # Concatenate forward hit tables

    cat $annotated_hits_collected > matches.dmnd.annot.tsv

    # Extract best forward hit pairs & their alignment length

    awk 'BEGIN{OFS="\t"}; {print \$1,\$4,\$9-\$8, "forward"}' matches.dmnd.annot.tsv | \
    uniq >> \
    mixed_hits.txt

    # Determine best single protein-locus pair from both the forward & reciprocal searches

    sort -k2,2 -k3,3nr mixed_hits.txt | \
    sort -u -k2,2 | \
    tee >(grep reciprocal-nr | cut -f1 | sort | uniq > extract_from_reciprocal-nr.txt) | \
    tee >(grep reciprocal-rvdb | cut -f1 | sort | uniq > extract_from_reciprocal-rvdb.txt) | \
    tee >(grep forward | cut -f1 | sort | uniq > extract_from_forward.txt) > \
    best_pairs_partial.txt

    # Add true genomic start for later conversion

    cut -f2 best_pairs_partial.txt | sed 's/[:-]/\t/g' > coords.txt
    paste best_pairs_partial.txt coords.txt > best_pairs.txt

    # Extract best proteins

    seqtk subseq reciprocal_nr_subset.fasta extract_from_reciprocal-nr.txt | \
    sed 's/ .*//' > \
    best_hits.fasta

    seqtk subseq reciprocal_rvdb_subset.fasta extract_from_reciprocal-rvdb.txt | \
    sed 's/ .*//' >> \
    best_hits.fasta

    seqtk subseq $clustered_proteins extract_from_forward.txt | \
    sed 's/ .*//' >> \
    best_hits.fasta

    """

}
