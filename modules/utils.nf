process CONCATENATE_PUBLISH_TABLES {
    input:
    path collected_files
    val table_name  // Receive the table name as a string

    output:
    path "${table_name}"
    publishDir "${params.outdir}/sql", mode: "copy", pattern: "${table_name}"

    """
    cat ${collected_files} > '${table_name}'
    """
}