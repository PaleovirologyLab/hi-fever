process attempt_genewise_improvement {

    input:

    path original_genewise
    path reciprocal_hits
    path reciprocal_fasta
    path reciprocal_db
    path strict_fastas_collected
    path context_fastas_collected
    path python_path

    """

    export WISECONFIGDIR="\$CONDA_PREFIX/share/wise2/wisecfg"
    mkdir wise_tmp

    # Generate table of original best proteins per locus

    cut -f5,8 $original_genewise | \
    sort -k1,1 > \
    wise_tmp/original_pairs

    # Generate table of reciprocal best proteins per locus

    sort -k1,1 -k12,12nr $reciprocal_hits | \
    sort -u -k1,1 | \
    cut -f1-2 > \
    wise_tmp/reciprocal_pairs

    # Find loci whose best hit has changed & get nuc + protein identifiers + genomic coords & file to intersect with context coords

    comm -13 wise_tmp/original_pairs wise_tmp/reciprocal_pairs | \
    tee >(cut -f1 > wise_tmp/nuc_headers) | \
    tee >(cut -f2 | sort | uniq > wise_tmp/protein_accessions) | \
    tee >(cut -f1 | sed 's/:/\t/; s/-/\t/' > wise_tmp/genomic_coords) | \
    tee >(sed 's/:/\t/; s/-/\t/' | sort -k1,1 -k2,2n > wise_tmp/intersection_bed) > \
    wise_tmp/matched_pairs

    # Initialise nextflow file path (due to nextflow bug?)

    touch $reciprocal_fasta

    # Extract new best proteins into individual files

    while read line
        do
            echo \$line | \
            seqtk subseq $reciprocal_fasta - > \
            wise_tmp/\$line
        done < wise_tmp/protein_accessions

    rm "\$(readlink -f $reciprocal_fasta)"

    # Get all strict FASTAs together and extract relevant ones

    cat $strict_fastas_collected | \
    sed 's/ .*//' > \
    wise_tmp/temp.fa

    while read line
        do
            echo \$line | \
            seqtk subseq wise_tmp/temp.fa - > \
            wise_tmp/\$line
        done < wise_tmp/nuc_headers

    # GeneWise operations, strict FASTA

    while read line
        do
            query=\$(echo \$line | cut -f2 -d ' ')
            target=\$(echo \$line | cut -f1 -d ' ')
            genewise wise_tmp/\$query wise_tmp/\$target -both -matrix "$params.genewise_matrix".bla -sum -pep -cdna -divide DIVIDE_STRING -silent | \
            grep -v ">\\|Bits   Query" | \
            awk '/^[-0-9]/ {printf("%s%s\\t",(N>0?"\\n":""),\$0);N++;next;} {printf("%s",\$0);} END {printf("\\n");}' | \
            sed 's/DIVIDE_STRING/\t/g' | \
            tr -s '\t' | \
            tr -s ' ' '\t' | \
            sort -k5,5 -k1,1nr | \
            sort -u -k5,5 >> \
            wise_tmp/genewise_strict
        done < wise_tmp/matched_pairs

    # Back-calculate genomic coords from genewise

    paste wise_tmp/genewise_strict wise_tmp/genomic_coords | \
    tr -s '\t' | \
    awk 'BEGIN{OFS="\t"} {if (\$6 < \$7) print \$12, \$13+\$6-1, \$13+\$7-1, "+", \$5, "strict_PR", \$1, \$2, \$3, \$4, \$8, \$9, \$10, \$11; else print \$12, \$13+\$7-1, \$13+\$6-1, "-", \$5, "strict_PR", \$1, \$2, \$3, \$4, \$8, \$9, \$10, \$11}' | \
    sort -k1,1 -k2,2n > \
    wise_tmp/output1

    # Get all context FASTAs together

    cat $context_fastas_collected | \
    sed 's/ .*//' > \
    wise_tmp/temp.fa

    #  Convert context ranges to BED for intersection with those strict regions where the best protein changed

    grep ">" wise_tmp/temp.fa | \
    sed 's/>//; s/:/\t/; s/-/\t/' | \
    sort -k1,1 -k2,2n > \
    wise_tmp/all_context_bed

    # Intersect ranges, keep various outputs from the context ranges we need

    bedtools intersect -a wise_tmp/intersection_bed -b wise_tmp/all_context_bed -wb -f 1 -sorted | \
    tee >(awk 'BEGIN{OFS="\t"} {print \$5":"\$6"-"\$7}' > wise_tmp/nuc_headers) | \
    tee >(cut -f5-7 > wise_tmp/genomic_coords) | \
    awk 'BEGIN{OFS="\t"} {print \$5":"\$6"-"\$7, \$4}' > \
    wise_tmp/matched_pairs

    # Extract those sequences

    while read line
        do
            echo \$line | \
            seqtk subseq wise_tmp/temp.fa - > \
            wise_tmp/\$line
        done < wise_tmp/nuc_headers

    # GeneWise operations, context FASTA

    while read line
        do
            query=\$(echo \$line | cut -f2 -d ' ')
            target=\$(echo \$line | cut -f1 -d ' ')
            genewise wise_tmp/\$query wise_tmp/\$target -both -matrix "$params.genewise_matrix".bla -sum -pep -cdna -divide DIVIDE_STRING -silent | \
            grep -v ">\\|Bits   Query" | \
            awk '/^[-0-9]/ {printf("%s%s\\t",(N>0?"\\n":""),\$0);N++;next;} {printf("%s",\$0);} END {printf("\\n");}' | \
            sed 's/DIVIDE_STRING/\t/g' | \
            tr -s '\t' | \
            tr -s ' ' '\t' | \
            sort -k5,5 -k1,1nr | \
            sort -u -k5,5 >> \
            wise_tmp/genewise_context
        done < wise_tmp/matched_pairs

    # Back-calculate genomic coords from genewise & remove redundancy (sites found in two context FASTAs)

    paste wise_tmp/genewise_context wise_tmp/genomic_coords | \
    tr -s '\t' | \
    awk 'BEGIN{OFS="\t"} {if (\$6 < \$7) print \$12, \$13+\$6-1, \$13+\$7-1, "+", \$5, "context_PR", \$1, \$2, \$3, \$4, \$8, \$9, \$10, \$11; else print \$12, \$13+\$7-1, \$13+\$6-1, "-", \$5, "context_PR", \$1, \$2, \$3, \$4, \$8, \$9, \$10, \$11}' | \
    sort -u -k1,1 -k2,2n -k3,3nr > \
    wise_tmp/output2

    # Intersect and concatenate results (keep strict if not covered by context, otherwise keep context)

    # First report strict predictions not encompassed by a context prediction
    bedtools intersect -v -a wise_tmp/output1 -b wise_tmp/output2 -f 1 -wa > wise_tmp/merged_results

    # Second find strict predictions encompassed by context predictions, report latter
    bedtools intersect -a wise_tmp/output1 -b wise_tmp/output2 -f 1 -wb | \
    awk 'BEGIN{OFS="\t"} {print \$15, \$16, \$17, \$18, \$5, \$20, \$21, \$22, \$23, \$24, \$25, \$26, \$27, \$28}' >> \
    wise_tmp/merged_results

    # Post-processing of in-frame STOPs

    python $python_path --task $params.stop_task --file wise_tmp/merged_results > post_reciprocal_genewise.txt

    # Merge both genewise outputs - for each locus keep the longest individual prediction

    cat pre_reciprocal_genewise.txt post_reciprocal_genewise.txt | \
    sort -k1,1 -k2,2n -k3,3nr | \
    sort -u -k5,5 > \
    genewise.txt

    # Cleanup

    rm -r wise_tmp
    rm "\$(readlink -f pre_reciprocal_genewise.txt)"  post_reciprocal_genewise.txt

    """

}