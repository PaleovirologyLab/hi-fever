process download_assemblies {

    maxForks 10

    input:
    path ftp_dir

    output:
    path "*genomic.fna.gz"

    """

    max_attempts=5

    count=0

    while read line
    do

    # Downloads and checks assembly file for corruption, re-attempts if md5 check fails
    md5check_function () {
        assemblyFile="\$(echo "\$line" | sed 's/^.*\\///')_genomic.fna.gz"
        wget -q "\$line/\$assemblyFile" "\$line/md5checksums.txt"
        grep \$assemblyFile md5checksums.txt > \$assemblyFile.md5; rm md5checksums.txt
        status=`md5sum -c \$assemblyFile.md5 2>/dev/null | sed 's/.* //'`
        if [ "\$status" == FAILED ]
        then
                        if [ "\$count" == "\$max_attempts" ]
                        then
                                echo "\$assemblyFile FAILED md5check \$max_attempts times, exiting"; exit 1
                        else
                                echo "\$assemblyFile FAILED md5check"; rm \$assemblyFile*; count=\$count +1; md5check_function
                        fi
        else
                echo "\$assemblyFile PASSED md5check"; rm \$assemblyFile.md5
        fi
        }

    md5check_function

    done < $ftp_dir

    """
}