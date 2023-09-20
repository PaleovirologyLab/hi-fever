process hmmer {

    input:
    path profile_dir
    path clustered_fasta

    output:
    path "query_domains.hmmer", emit: query_domains_ch
    publishDir "query_domains"

    """

    hmmscan --noali --notextw --qformat fasta --domtblout raw_domains.txt $profile_dir/*.hmm $clustered_fasta 1> /dev/null

    # Post-processing:
    # Merge overlapping query protein alignments, keep best (by bitscore)
    # For each best alignment, if e-value is under 0.1, report:
    # query overlap_start overlap_end best_start best_end best_bitscore best_i-Evalue model_acc model_name model_description

    grep -v "#" raw_domains.txt | \
    tr -s ' ' '\t' | \
    awk 'BEGIN{OFS="\t"}; {print \$4, \$18, \$19, \$13, \$14, \$2, \$1, \$23}' | \
    sort -k1,1 -k2,2n | \
    tee >(bedtools merge > tmp.out) | \
    bedtools cluster | \
    sort -k9,9n -k5,5nr | \
    sort -u -k9,9n | \
    bedtools intersect -a - -b tmp.out -wb -sorted | \
    awk 'BEGIN{OFS="\t"}; {if (\$4 < 0.1) print \$1, \$11, \$12, \$2, \$3, \$5, \$4, \$6, \$7, \$8}' > \
    query_domains.hmmer

    rm tmp.out

    """

}