process diamond {

    maxForks params.diamond_forks

    cpus 12
    machineType 'n1*'
    memory '200 GB'
    disk '100 GB'

    input:
    tuple path(assembly), path(db), val(ready)

    output:
    tuple path("*.dmnd.tsv"), path("*.gz*nsq")

    """

    db=$db
    assembly=$assembly

    diamond blastx -v \
    --$params.diamond_mode \
    --matrix $params.diamond_matrix \
    --masking seg \
    -d \$db \
    -q \$assembly \
    -o matches.dmnd.tsv \
    --outfmt 6 qseqid qstart qend qframe qlen sseqid sstart send slen evalue bitscore pident length mismatch gapopen &

    gunzip -c \$assembly | makeblastdb -in - -out \${assembly/*\\//} -title \${assembly/*\\//} -dbtype nucl -parse_seqids &

    wait

    rm "\$(readlink -f \$assembly)"

    """

}
