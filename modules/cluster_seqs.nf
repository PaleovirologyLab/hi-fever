process CLUSTER_SEQS {

	container 'oras://community.wave.seqera.io/library/mmseqs2:17.b804f--7248ed4d147f7327'
	conda 'bioconda::mmseqs2=17.b804f'

    input:
    path queries

    output:
    path "query_db_clu_rep.fasta"
    publishDir "${params.outdir}", mode: "copy"

    """

    mmseqs createdb $queries query_db
    mmseqs linclust --min-seq-id $params.mmseqs_minseqid --cov-mode 1 -c $params.mmseqs_cover query_db query_db_clu tmp
    mmseqs createsubdb query_db_clu query_db query_db_clu_rep
    mmseqs convert2fasta query_db_clu_rep query_db_clu_rep.fasta
    rm -r tmp

    """

}
