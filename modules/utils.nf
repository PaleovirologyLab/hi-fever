process CHECK_FILE_TYPE {
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

process CONCATENATE_PUBLISH_TABLES {
    input:
    path collected_files
    val table_name  // Receive the table name as a string

    output:
    path "${table_name}"
    publishDir "${params.outdir}/sql", mode: "copy", pattern: "${table_name}"

    """
    cat $collected_files > '${table_name}'
    """
}