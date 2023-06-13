#!/home/user/nextflow

// Syntax version
nextflow.enable.dsl=2

// Script parameters
params.ftp_file = "$PWD/ftp_list.txt"
params.query_file_aa = "$PWD/protein_query.fasta"
params.mmseqs_minseqid = "0.95"
params.mmseqs_cover = "0.90"
params.diamond_mode = "very-sensitive"
params.diamond_matrix = "BLOSUM62"
params.diamond_cpus = "12"

process build_db {

    input:
    path x

    output:
    path "DB_clu_rep.fasta"
    path "virusdb.dmnd"
    publishDir "virusdb"

    """

    mmseqs createdb $x DB
    mmseqs cluster --min-seq-id $params.mmseqs_minseqid -c $params.mmseqs_cover DB DB_clu tmp
    mmseqs createsubdb DB_clu DB DB_clu_rep
    mmseqs convert2fasta DB_clu_rep DB_clu_rep.fasta
    diamond makedb --in DB_clu_rep.fasta -d virusdb

    """

}

process parse_ftp {

    input:
    path x

    output:
    path "*.ftp.txt"
    

    """

    while read line
    do

    assembly_ID=`echo \$line | sed 's/^.*\\///'`

    echo \$line > \${assembly_ID}.ftp.txt


    done < $x

    """


}

process download_assemblies {

    //debug true

    input:
    path y

    output:
    path "*genomic.fna.gz"

    """

    max_downloads=5

    count=0

    while read line
    do

    # Downloads and checks assembly file for corruption, re-attempts if md5 check fails
    md5check_function () {
        assembly_ID=`echo \$line | sed 's/^.*\\///'`
        wget -q \$line/\${assembly_ID}_genomic.fna.gz
        wget -q \$line/md5checksums.txt
        assemblyFile=`ls -1 *genomic.fna.gz`
        grep \$assemblyFile md5checksums.txt > \$assemblyFile.md5; rm md5checksums.txt
        status=`md5sum -c \$assemblyFile.md5 2>/dev/null | sed 's/.* //'`
        if [ "\$status" == FAILED ]
        then
                        if [ "\$count" == "\$max_downloads" ]
                        then
                                echo "\$assemblyFile FAILED md5check \$max_downloads times, exiting"; exit 1
                        else
                                echo "\$assemblyFile FAILED md5check"; rm \$assemblyFile*; count=\$count +1; md5check_function
                        fi
        else
                echo "\$assemblyFile PASSED md5check"; rm \$assemblyFile.md5
        fi
        }

    md5check_function

    done < $y

    """
}

process diamond {

    maxForks 4
    //debug true

    input:
    path x

    output:
    path "*.dmnd.tsv"

    """

    diamond blastx \
    --$params.diamond_mode \
    --matrix $params.diamond_matrix \
    --masking seg \
    -d $PWD/virusdb/virusdb.dmnd \
    -q $x \
    -o matches.dmnd.tsv \
    -p $params.diamond_cpus \
    --outfmt 6 qseqid qstart qend qframe qlen sseqid sstart send slen evalue bitscore pident length mismatch gapopen

    """

}

process bedtools {

    input:
    path x

    output:
    path "*.nonredundant.bed"

    """

    awk 'BEGIN{OFS="\t"}; {if(\$2<\$3) print \$1, \$2, \$3; else if(\$3<\$2) print \$1, \$3, \$2}' \
    $x | sort -k1,1 -k2,2n | bedtools merge -c 1 -o count > diamond-result.nonredundant.bed

    """

}

workflow {
    def db_ch = Channel.fromPath(params.query_file_aa)
    build_db(db_ch)
    def ftp_ch = Channel.fromPath(params.ftp_file)
    parse_ftp(ftp_ch).flatten() | download_assemblies | diamond | bedtools

}
