process diamond {

    maxForks params.diamond_forks

    input:
    tuple path(assembly), path(db)

    output:
    tuple path("*.dmnd.tsv"), path("*.gz*nsq")

    """

    db=$db
    assembly=$assembly
    chunks=\$(echo \$assembly | sed 's/_genomic.*/_genomic_chunks.fna.gz/')
    cpu_count=\$(awk -v total_cpu=\$(nproc) 'BEGIN {printf "%.0f\\n", (total_cpu > 1) ? total_cpu / $params.diamond_forks : 1}')

    seqkit sliding -s $params.chunk_size -W $params.chunk_size -g $assembly -o \$chunks

    diamond blastx \
    --$params.diamond_mode \
    --matrix $params.diamond_matrix \
    --masking seg \
    -d \$db \
    -q \$chunks \
    -o matches.out \
    -p \$cpu_count \
    --max-target-seqs 1000 \
    --outfmt 6 qseqid qstart qend qframe qlen sseqid sstart send slen evalue bitscore pident length mismatch gapopen &

    gunzip -c \$assembly | makeblastdb -in - -out \${assembly/*\\//} -title \${assembly/*\\//} -dbtype nucl -parse_seqids &

    wait

    sed 's/_sliding:/\\t/' matches.out | sed 's/-/\\t/' | awk -v OFS='\\t' '{print \$1,\$2+\$4-1,\$2+\$5-1,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15,\$16,\$17}' > matches.reformatted.dmnd.tsv

    rm \$chunks

    rm "\$(readlink -f \$assembly)"

    """

}