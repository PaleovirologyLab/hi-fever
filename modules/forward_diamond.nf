process FORWARD_DIAMOND {

    maxForks params.diamond_forks
	tag "${meta.id}"

	container 'oras://community.wave.seqera.io/library/diamond_seqkit_seqtk:6fc81cc10da8e7e4'
	conda 'bioconda::diamond=2.1.11 bioconda::seqkit=2.9.0'

    input:
    tuple val(meta), path(assembly), path(db)

    output:
    tuple val(meta), path("*.dmnd.tsv")

    //publishDir "${params.outdir}/forwardDiamond", mode: "copy", pattern: "*.dmnd.tsv"

    """

    chunks=\$(echo ${assembly} | sed 's/_genomic.*/_genomic_chunks.fna.gz/')
    cpu_count=\$(awk -v total_cpu=\$(nproc) 'BEGIN {printf "%.0f\\n", (total_cpu > 1) ? total_cpu / ${params.diamond_forks} : 1}')

    seqkit sliding -s ${params.chunk_size} -W ${params.chunk_size} -g ${assembly} -o \$chunks

    diamond blastx \
    --${params.diamond_mode} \
    --matrix ${params.diamond_matrix} \
    --masking seg \
    -d ${db} \
    -q \$chunks \
    -o matches.out \
    -p \$cpu_count \
    --max-target-seqs ${params.diamond_max_target_seqs} \
    --outfmt 6 qseqid qstart qend qframe qlen sseqid sstart send slen evalue bitscore pident length mismatch gapopen

	sed 's/_sliding:/\\t/' matches.out | sed 's/-/\\t/' | \
		awk -v OFS='\\t' '{print \$1,\$2+\$4-1,\$2+\$5-1,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15,\$16,\$17}' \
		> "${assembly}_forward-matches-raw.dmnd.tsv"

    """

}
