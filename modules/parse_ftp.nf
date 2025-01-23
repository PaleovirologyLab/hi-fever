process parse_ftp {

    input:
    path ftp_input

    output:
    path "*.ftp.txt"

    """

    while read line
        do
            assemblyID=`echo \$line | sed 's/^.*\\///'`
            echo \$line > \${assemblyID}.ftp.txt
        done < $ftp_input

    """

}

process check_file_type {
    input:
    path inputFile // Input file

    output:
    tuple path(inputFile), val(type) // Output the file and its type (fasta or dmnd)

    script:
    def extension = inputFile.name.tokenize('.')[-1]

    if (extension == 'fasta') {
        type = 'fasta'
    } else if (extension == 'fa') {
        type = 'fasta'
    } else if (extension == 'dmnd') {
        type = 'dmnd'
    } else {
        error "Unsupported file extension: ${extension}"
    }

    """
    echo "$type"
    """
}