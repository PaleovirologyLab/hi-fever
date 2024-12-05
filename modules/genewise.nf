process genewise {

    input:
    path pair_subsets
    path best_hit_proteins
    path strict_fastas
    path context_fastas
    path context_coords_val

    output:
    stdout

    """

    # General setup for genewise

    export WISECONFIGDIR="\$CONDA_PREFIX/share/wise2/wisecfg"
    mkdir wise_tmp

    # Prepare lists with sequences required per task, their pairing (loop input), original genomic coordinates, & BED of coords to intersect with context coords

    cut -f1 $pair_subsets | sort | uniq > wise_tmp/required_proteins.txt
    cut -f2 $pair_subsets | sort | uniq > wise_tmp/required_strict_nts.txt
    cut -f1-2 $pair_subsets > wise_tmp/loop_input1.txt
    cut -f5-7 $pair_subsets > wise_tmp/genomic_coords1
    awk 'BEGIN {OFS="\t"} {print \$5, \$6, \$7, \$1}' $pair_subsets | sort -k1,1 -k2,2n > wise_tmp/intersection_bed

    # Extract required sequence subsets

    seqtk subseq $best_hit_proteins wise_tmp/required_proteins.txt | \
    sed 's/ .*//' > \
    wise_tmp/subset_proteins.fa

    seqtk subseq $strict_fastas wise_tmp/required_strict_nts.txt | \
    sed 's/ .*//' > \
    wise_tmp/subset_strict_nts.fa

    # Process subsets into individual files

    awk '/^>/ { if (name) close(name); name="wise_tmp/" substr(\$0,2); print > name; next } { print >> name }' wise_tmp/subset_proteins.fa

    awk '/^>/ { if (name) close(name); name="wise_tmp/" substr(\$0,2); print > name; next } { print >> name }' wise_tmp/subset_strict_nts.fa

    # GeneWise operations, strict FASTA

    while read line
        do
            query=\$(echo \$line | cut -f1 -d ' ')
            target=\$(echo \$line | cut -f2 -d ' ')
            genewise wise_tmp/\$query wise_tmp/\$target -both -kbyte 2000000 -matrix "$params.genewise_matrix".bla -sum -cdna -divide DIVIDE_STRING -silent | \
            grep -v ">\\|Bits   Query\\|intron" | \
            awk '/^[-0-9]/ {printf("%s%s\\t",(N>0?"\\n":""),\$0);N++;next;} {printf("%s",\$0);} END {printf("\\n");}' | \
            sed 's/DIVIDE_STRING/\t/g' | \
            tr -s '\t' | \
            tr -s ' ' '\t' | \
            sort -k5,5 -k1,1nr | \
            sort -u -k5,5 >> \
            wise_tmp/genewise_strict
        done < wise_tmp/loop_input1.txt

    # Back-calculate genomic coords from genewise

    paste wise_tmp/genewise_strict wise_tmp/genomic_coords1 | \
    tr -s '\t' | \
    awk 'BEGIN{OFS="\t"} {if (\$6 < \$7) print \$11, \$12+\$6-1, \$12+\$7-1, "+", \$5, "strict", \$1, \$2, \$3, \$4, \$10, \$9, \$8; else print \$11, \$12+\$7-1, \$12+\$6-1, "-", \$5, "strict", \$1, \$2, \$3, \$4, \$10, \$9, \$8}' | \
    sort -k1,1 -k2,2n > \
    wise_tmp/output1

    # Find any context regions fully overlapping strict regions in this task

    bedtools intersect -a wise_tmp/intersection_bed -b all_context_coords.bed -loj -f 1 -sorted | \
    tee wise_tmp/genomic_coords2 | \
    awk 'BEGIN {OFS="\t"} {print \$4, \$5":"\$6"-"\$7}' > \
    wise_tmp/loop_input2.txt

    # Prepare list of context sequences required for this task

    cut -f2 wise_tmp/loop_input2.txt | sort | uniq > wise_tmp/required_context_nts.txt

    # Extract required sequence subset

    seqtk subseq $context_fastas wise_tmp/required_context_nts.txt | \
    sed 's/ .*//' > \
    wise_tmp/subset_context_nts.fa

    # Process subsets into individual files

    awk '/^>/ { if (name) close(name); name="wise_tmp/" substr(\$0,2); print > name; next } { print >> name }' wise_tmp/subset_context_nts.fa

    # GeneWise operations, context FASTA

    while read line
        do
            query=\$(echo \$line | cut -f1 -d ' ')
            target=\$(echo \$line | cut -f2 -d ' ')
            genewise wise_tmp/\$query wise_tmp/\$target -both -kbyte 2000000 -matrix "$params.genewise_matrix".bla -sum -cdna -divide DIVIDE_STRING -silent | \
            grep -v ">\\|Bits   Query\\|intron" | \
            awk '/^[-0-9]/ {printf("%s%s\\t",(N>0?"\\n":""),\$0);N++;next;} {printf("%s",\$0);} END {printf("\\n");}' | \
            sed 's/DIVIDE_STRING/\t/g' | \
            tr -s '\t' | \
            tr -s ' ' '\t' | \
            sort -k5,5 -k1,1nr | \
            sort -u -k5,5 >> \
            wise_tmp/genewise_context
        done < wise_tmp/loop_input2.txt

    # Back-calculate genomic coords from genewise, remove redundancy (sites found in two context FASTAs)

    paste wise_tmp/genewise_context wise_tmp/genomic_coords2 | \
    tr -s '\t' | \
    awk 'BEGIN{OFS="\t"} {if (\$6 < \$7) print \$15, \$16+\$6-1, \$16+\$7-1, "+", \$11":"\$12"-"\$13, "context", \$1, \$2, \$3, \$4, \$10, \$9, \$8; else print \$15, \$16+\$7-1, \$16+\$6-1, "-", \$11":"\$12"-"\$13, "context", \$1, \$2, \$3, \$4, \$10, \$9, \$8}' | \
    sort -u -k5,5 | \
    sort -k1,1 -k2,2n > \
    wise_tmp/output2

    # Report strict regions not encompassed by a context prediction

    bedtools intersect -v -a wise_tmp/output1 -b wise_tmp/output2 -f 1 -wa | \
    cut -f5 > \
    wise_tmp/remove_from_context

    # Remove off-target context alignments (not covering the target locus)
    # N.B. an empty grep will produce non-zero exit status which Nextflow abhors
    # "|| true" forces a 0 exit if that occurs

    grep -vFf wise_tmp/remove_from_context wise_tmp/output2 > wise_tmp/output2_on_target || true

    # Merge results

    cat wise_tmp/output1 wise_tmp/output2_on_target | \
    sort -k5,5 -k2,2n -k3,3nr | \
    sort -u -k5,5 > \
    wise_tmp/merged_results
    """

}

    // # Post-processing of in-frame STOPs

    // stopConvertAndCount.py --task $params.stop_task --file wise_tmp/merged_results > wise_tmp/stops_processed

    // # Translation of coding sequences
    // # Outputs: contig genomic_start genomic_end strand locus sourceFASTA bitscore query qstart qend cdna peptide intron_count idels_frameshifts inframe_STOPs

    // translateCodingSequence.py --input wise_tmp/stops_processed --output /dev/stdout
