process genewise {

    input:
    tuple path(annotated_tsv), path(strict_fasta), path(context_fasta), path(context_coords), path(clustered_proteins), path(python_path)

    output:
    path "*_genewise"

    """

    # Check if any candidates were found in the assembly to proceed

    title=\$(readlink -f $strict_fasta | sed 's/.*\\///; s/_genomic.fna_strict.fasta//')

    if [ -s $strict_fasta ];

        then

            protein_db=$clustered_proteins
            export WISECONFIGDIR="\$CONDA_PREFIX/share/wise2/wisecfg"
            mkdir wise_tmp

            # Reformat and extract strict nucleotide FASTA, one per file

            sed 's/ .*//' < $strict_fasta > wise_tmp/temp.fa
            grep ">" wise_tmp/temp.fa | sed 's/>//' > wise_tmp/nuc_headers

            while read line
                do
                    echo \$line | \
                    seqtk subseq wise_tmp/temp.fa - > \
                    wise_tmp/\$line
                done < wise_tmp/nuc_headers

            # Generate query-target pairing file for strict FASTA, protein accessions to extract, file for later intersection with context coordinates, and genomic coordinates file
            # First cut protein, contig, strict coords start and end -> send to matched_pairs

            cut -f1,4-6 $annotated_tsv | \
            uniq | \
            tee >(cut -f1 | sort | uniq > wise_tmp/prot_headers) | \
            tee >(awk 'BEGIN{OFS="\t"} {print \$2, \$3, \$4, \$1}' | sort -k1,1 -k2,2n > wise_tmp/intersection_bed) | \
            tee >(cut -f2-4 > wise_tmp/genomic_coords) | \
            awk 'BEGIN{OFS="\t"} {print \$1,\$2":"\$3"-"\$4}' > \
            wise_tmp/matched_pairs

            # Extract protein FASTA, one per file

            while read line
                do
                    echo \$line | \
                    seqtk subseq \$protein_db - > \
                    wise_tmp/\$line
                done < wise_tmp/prot_headers

            # GeneWise operations, strict FASTA

            while read line
                do
                    query=\$(echo \$line | cut -f1 -d ' ')
                    target=\$(echo \$line | cut -f2 -d ' ')
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
            awk 'BEGIN{OFS="\t"} {if (\$6 < \$7) print \$12, \$13+\$6-1, \$13+\$7-1, "+", \$5, "strict", \$1, \$2, \$3, \$4, \$8, \$9, \$10, \$11; else print \$12, \$13+\$7-1, \$13+\$6-1, "-", \$5, "strict", \$1, \$2, \$3, \$4, \$8, \$9, \$10, \$11}' | \
            sort -k1,1 -k2,2n > \
            wise_tmp/output1

            # Reformat and extract context nucleotide FASTA, one per file

            sed 's/ .*//' < $context_fasta > wise_tmp/temp.fa
            grep ">" wise_tmp/temp.fa | sed 's/>//' > wise_tmp/nuc_headers

            while read line
                do
                    echo \$line | \
                    seqtk subseq wise_tmp/temp.fa - > \
                    wise_tmp/\$line
                done < wise_tmp/nuc_headers

            # Generate query-target pairing file for context FASTA

            bedtools intersect -a $context_coords -b wise_tmp/intersection_bed -loj -F 1 -sorted | \
            cut -f1-3,7 | \
            tee wise_tmp/genomic_coords | \
            awk 'BEGIN{OFS="\t"} {print \$4, \$1":"\$2"-"\$3}' > \
            wise_tmp/matched_pairs

            # GeneWise operations, context FASTA

            while read line
                do
                    query=\$(echo \$line | cut -f1 -d ' ')
                    target=\$(echo \$line | cut -f2 -d ' ')
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

            # Back-calculate genomic coords from genewise, remove redundancy (sites found in two context FASTAs)

            paste wise_tmp/genewise_context wise_tmp/genomic_coords | \
            tr -s '\t' | \
            awk 'BEGIN{OFS="\t"} {if (\$6 < \$7) print \$12, \$13+\$6-1, \$13+\$7-1, "+", \$5, "context", \$1, \$2, \$3, \$4, \$8, \$9, \$10, \$11; else print \$12, \$13+\$7-1, \$13+\$6-1, "-", \$5, "context", \$1, \$2, \$3, \$4, \$8, \$9, \$10, \$11}' | \
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

            python $python_path --task $params.stop_task --file wise_tmp/merged_results > "\${title}_genewise"

            # Cleanup

            rm -r wise_tmp

        else

            touch "\${title}_genewise"

    fi

    """

}
