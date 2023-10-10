process build_db {

    input:
    path queries

    output:
    path "DB_clu_rep.fasta", emit: clust_ch
    path "virusdb.dmnd", emit: vir_db_ch
    publishDir "${params.outdir}/virusdb", mode: "copy"


    """

    mmseqs createdb $queries DB
    mmseqs cluster --min-seq-id $params.mmseqs_minseqid --cov-mode 1 -c $params.mmseqs_cover DB DB_clu tmp
    mmseqs createsubdb DB_clu DB DB_clu_rep
    mmseqs convert2fasta DB_clu_rep DB_clu_rep.fasta
    diamond makedb --in DB_clu_rep.fasta -d virusdb
    rm -r tmp

    """

}