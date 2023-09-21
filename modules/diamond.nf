process diamond {

    maxForks params.diamond_forks

    input:
    tuple path(assembly), path(db)

    output:
    tuple path("*.dmnd.tsv"), path("*.gz*nsq")

    """

    db=$db
    assembly=$assembly
    cpu_count=\$(awk -v total_cpu=\$(nproc) 'BEGIN {printf "%.0f\\n", (total_cpu > 1) ? total_cpu / $params.diamond_forks : 1}')

    diamond blastx \
    --$params.diamond_mode \
    --matrix $params.diamond_matrix \
    --masking seg \
    -d \$db \
    -q \$assembly \
    -o matches.dmnd.tsv \
    -p \$cpu_count \
    --outfmt 6 qseqid qstart qend qframe qlen sseqid sstart send slen evalue bitscore pident length mismatch gapopen &

    gunzip -c \$assembly | makeblastdb -in - -out \${assembly/*\\//} -title \${assembly/*\\//} -dbtype nucl -parse_seqids &

    wait

    rm "\$(readlink -f \$assembly)"

    """

}