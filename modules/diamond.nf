process diamond {

    maxForks 4

    input:
    tuple path(assembly), path(db)

    output:
    tuple path("*.dmnd.tsv"), path("*.gz*nsq")

    """

    db=$db
    assembly=$assembly

    diamond blastx \
    --$params.diamond_mode \
    --matrix $params.diamond_matrix \
    --masking seg \
    -d \$db \
    -q \$assembly \
    -o matches.dmnd.tsv \
    -p $params.diamond_cpus \
    --outfmt 6 qseqid qstart qend qframe qlen sseqid sstart send slen evalue bitscore pident length mismatch gapopen &

    gunzip -c \$assembly | makeblastdb -in - -out \${assembly/*\\//} -title \${assembly/*\\//} -dbtype nucl -parse_seqids &

    wait

    rm "\$(readlink -f \$assembly)"

    """

}