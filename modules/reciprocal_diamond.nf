process reciprocal_diamond {

    input:
    path genewise
    path reciprocal_db

    output:
    path "reciprocal-matches.dmnd.tsv", emit: reciprocal_hits_ch
    path "pre_reciprocal_genewise.txt", emit: original_genewise_ch
    path "best_reciprocals.fasta", emit: reciprocal_fasta_ch
    publishDir "${params.outdir}/sql", mode: "copy", pattern: "reciprocal-matches.dmnd.tsv"

    """

    cat $genewise > pre_reciprocal_genewise.txt

    awk '{print ">"\$5"\\n"\$13}' pre_reciprocal_genewise.txt > predicted_proteins.fasta

    diamond blastp \
    --$params.diamond_mode \
    --matrix $params.diamond_matrix \
    --masking seg \
    -d $reciprocal_db \
    -q predicted_proteins.fasta \
    -e 1e-5 \
    -k 20 \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms skingdoms sphylums full_sseq | \
    tee >(sort -k1,1 -k12,12nr | sort -u -k1,1 | cut -f2,18 | sort -u -k1,1 | awk '{print ">"\$1"\\n"\$2}' > best_reciprocals.fasta) | \
    cut -f1-17 > \
    reciprocal-matches.dmnd.tsv

    """

}