process hmmer {

    input:
    path profile_dir
    path clustered_fasta

    output:
    path "query_domains.hmmer", emit: query_domains_ch
    publishDir "query_domains"

    """

    mkdir tmp

    # Get domain descriptions to add later

    grep "ACC \\|DESC " $profile_dir/*.hmm | \
    tr -s ' ' '\t' | \
    cut -f2 | \
    awk 'NR%2{printf "%s\t",\$0;next;}1' > \
    tmp/descriptions

    # Run HMM search

    hmmsearch --noali --notextw --domtblout tmp/raw_domains.txt $profile_dir/*.hmm $clustered_fasta 1> /dev/null

    # Merge overlapping query protein alignments, keep best (by bitscore)

    grep -v "#" tmp/raw_domains.txt | \
    tr -s ' ' '\t' | \
    awk 'BEGIN{OFS="\t"}; {print \$1, \$18, \$19, \$13, \$14, \$5, \$4}' | \
    sort -k1,1 -k2,2n | \
    tee >(bedtools merge > tmp/tmp1.out) | \
    bedtools cluster | \
    sort -k8,8n -k5,5nr | \
    sort -u -k8,8n | \
    tee >(cut -f6 > tmp/accessions) | \
    bedtools intersect -a - -b tmp/tmp1.out -wb -sorted > \
    tmp/tmp2.out

    # For each domain accession found, get description

    while read line
        do
            grep \$line tmp/descriptions | \
            cut -f2 >> \
            tmp/desc_merge
        done < tmp/accessions

    # Merge in domain description
    # For each best alignment, if e-value is under 0.1, report:
    # query overlap_start overlap_end best_start best_end best_bitscore best_i-Evalue model_acc model_name model_description

    paste tmp/tmp2.out tmp/desc_merge | \
    awk 'BEGIN{OFS="\t"}; {if (\$4 < 0.1) print \$1, \$10, \$11, \$2, \$3, \$5, \$4, \$6, \$7, \$12}' > \
    query_domains.hmmer

    rm -r tmp

    """

}