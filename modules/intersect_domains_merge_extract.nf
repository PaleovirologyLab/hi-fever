process extract_seqs_annotate_matches {

    input:
    tuple path(diamond_tsv), path(assembly_nsq_db)

    output:
    path "strict_coords.bed", emit: strict_coords_ch
    path "*_matches.dmnd.annot.tsv", emit: annot_tsv_ch
    path "*_strict.fasta", emit: strict_fa_ch
    path "*_context.fasta", emit: context_fa_ch
    path "*_locus_assembly_map.tsv", emit: locus_assembly_map_ch

    """
    # Prepare variables

    dbpath=\$(readlink -f \$(echo $assembly_nsq_db | cut -d ' ' -f1) | sed 's/.nsq//g; s/\\.[0-9][0-9]\$//g')
    filename=\$(echo \$dbpath | sed 's/\\.gz//g; s/\\/.*\\///g')
    assemblyID=\$(echo \$filename | sed 's/_genomic.*//')

    # Intersect domains & produce non-redundant BED:
    # Converts DIAMOND tsv to ascending assembly coordinate ranges, sorts to BED compatibility (contig and start position).
    # Calculates maximal strictly overlapping coordinate ranges and stores in BED file.
    # Adds a column to main tsv denoting which coordinates are in mergable overlapping clusters
    # Sorts best bitscore alignment per cluster to the top and removes others
    # Intersects this with the maximal range temp file
    # Converts to a subject oriented bed (i.e., protein hit and coordinates).

    # Reports:
    # sseqid_(protein) sstart send locus_id qseqid_best qstart_overlap qend_overlap qstart_best qend_best qframe_best qlen_best slen
    # eval bitscore pident len mismatch gapopen domain_overlap_start domain_overlap_end best_start best_end best_bitscore best_i-Evalue
    # model_acc model_name model_description

    awk 'BEGIN{OFS="\t"}; {if(\$2<\$3) print \$0; \
    else if (\$3<\$2) print \$1, \$3, \$2, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13, \$14, \$15}' $diamond_tsv | \
    sort -k1,1 -k2,2n | \
    tee >(bedtools merge > strict_coords.bed) | \
    bedtools cluster | \
    sort -k16,16n -k11,11nr | \
    sort -u -k16,16n | \
    bedtools intersect -a - -b strict_coords.bed -loj -sorted | \
    awk 'BEGIN{OFS="\t"}; \
    {print \$6, \$7, \$8, \$17":"\$18"-"\$19, \$17, \$18, \$19, \$2, \$3, \$4, \$5, \$9, \$10, \$11, \$12, \$13, \$14, \$15}' | \
    sort -k1,1 -k2,2n | \
    cut -f19 --complement > \
    \${assemblyID}_matches.dmnd.annot.tsv

    # First coordinate range extraction (strictly overlapping alignments)

    awk '{print \$1, \$2"-"\$3}' strict_coords.bed | \
    blastdbcmd -entry_batch - -db \$dbpath > \
    "\${filename}_strict.fasta"

    # Second coordinate range extraction (allow interval and add flanks)

    bedtools merge -d $params.interval -i strict_coords.bed | \
    awk -v flank=$params.flank '{if(\$2-flank < 1) print \$1, 1"-"\$3+flank; else print \$1, \$2-flank"-"\$3+flank}' | \
    blastdbcmd -entry_batch - -db \$dbpath > \
    "\${filename}_context.fasta"

    # Generate assemblyID to locus dictionary

    awk -v var="\$assemblyID" 'BEGIN{OFS="\t"}; {print \$1":"\$2"-"\$3, var}' strict_coords.bed > \
    "\${assemblyID}_locus_assembly_map.tsv"

    # Clean up database files

    rm \$dbpath*

    """

}