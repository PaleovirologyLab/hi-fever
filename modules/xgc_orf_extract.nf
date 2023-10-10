process orf_extract {

    input:
    path context_fasta
    path strict_coords

    output:
    path "context_ORFs.txt", optional: true
    path "intersected_ORFs.txt", optional: true

    """

    # Checks input file for content, if yes:
    # Extracts ORFs between START and STOP codons. STOP is NOT included in reported ORF coordinates.
    # Converts ORF output to the original genomic loci in ascending BED format, with predicted protein sequence.

    if [ -s $context_fasta ];

        then

            cat $context_fasta | \
            sed 's/:/-/' | \
            getorf -sequence /dev/stdin -outseq /dev/stdout -minsize $params.orf_size_nt -find 1 2>/dev/null | \
            seqtk seq - | \
            sed -r 's/] .*//; s/_[0-9]+ \\[/ /; s/>//; s/-/ /g' | \
            tr -s ' ' '\t' | \
            paste -sd '\t\n' | \
            awk 'BEGIN{OFS="\t"}; {if (\$4 < \$5) print \$1, \$2+\$4-1, \$2+\$5-1, "+", \$6; else print \$1, \$2+\$5-1, \$2+\$4-1, "-", \$6}' | \
            sort -k1,1 -k2,2n > \
            context_ORFs.txt

            # Intersect ORFs with strictly overlapping features, reports:
            # contig strict_feature_start strict_feature_end coverage_of_feature_by_ORF coverage_of_ORF_by_feature ORF_start ORF_end ORF_strand ORF_seq

            bedtools intersect -a $strict_coords -b context_ORFs.txt -wo -sorted | \
            awk '{print \$1, \$2, \$3, \$9/(\$3-\$2), (\$3-\$2)/(\$6-\$5), \$5, \$6, \$7, \$8}' | \
            tr -s ' ' '\t' > \
            intersected_ORFs.txt

        else

            true

    fi

    """

}