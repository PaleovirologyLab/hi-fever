nextflow.enable.dsl=2

include { CONCATENATE_PUBLISH_TABLES } from '../../modules/utils.nf'

workflow {
    files = Channel.fromPath(params.input_glob, checkIfExists: true).collect()
    CONCATENATE_PUBLISH_TABLES(files, params.table_name)
}
