nextflow.enable.dsl=2

include { CREATE_SUMMARY_TABLE_FULL } from '../../modules/create_summary_table.nf'

workflow {
    nr_ch = Channel.fromPath(params.reciprocal_nr, checkIfExists: true)
    rvdb_ch = Channel.fromPath(params.reciprocal_rvdb, checkIfExists: true)
    tax_ch = Channel.fromPath(params.taxonomy, checkIfExists: true)
    map_ch = Channel.fromPath(params.assembly_map, checkIfExists: true)
    meta_ch = Channel.fromPath(params.assembly_metadata, checkIfExists: true)
    gw_ch = Channel.fromPath(params.genewise, checkIfExists: true)

    out = CREATE_SUMMARY_TABLE_FULL(nr_ch, rvdb_ch, tax_ch, map_ch, meta_ch, gw_ch)

    out.summary_table
        .map { it.getName() }
        .collectFile(name: 'summary_manifest.txt', newLine: true, storeDir: params.outdir)

    out.aa_fasta
        .map { it.getName() }
        .collectFile(name: 'aa_manifest.txt', newLine: true, storeDir: params.outdir)

    out.cdna_fasta
        .map { it.getName() }
        .collectFile(name: 'cdna_manifest.txt', newLine: true, storeDir: params.outdir)
}
