#!/home/user/nextflow

// Syntax version
nextflow.enable.dsl=2

// Script parameters
params.ftp_file = "$PWD/ftp_list.txt"
params.query_file_aa = "$PWD/protein_query.fasta"
params.diamond_mode ="--very-sensitive"
params.diamond_matrix="BLOSUM62"

process build_diamond_db() {

    output:
    path "*.dmnd"
    publishDir "virusdb"

    """
    diamond makedb --in $params.query_file_aa -d virusdb

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

    debug true

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

    maxForks 1
    debug true

    input:
    path x

    output:
    path "*.dmnd.tsv"

    """
    
    diamond blastx $params.diamond_mode --matrix $params.diamond_matrix --masking seg -d $PWD/virusdb/virusdb.dmnd -q $x -o matches.dmnd.tsv

    """

}

workflow {
    build_diamond_db()
    def ftp_ch = Channel.fromPath(params.ftp_file)
    parse_ftp(ftp_ch).flatten() | download_assemblies | diamond
}