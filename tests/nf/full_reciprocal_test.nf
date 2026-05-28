nextflow.enable.dsl=2

include { FULL_RECIPROCAL_DIAMOND } from '../../modules/diamond.nf'

workflow {
    nr_db_ch = Channel.fromPath(params.reciprocal_nr_db, checkIfExists: true)
    rvdb_db_ch = Channel.fromPath(params.reciprocal_rvdb_db, checkIfExists: true)
    loci_ch = Channel.fromPath(params.loci_merged_fa, checkIfExists: true)
    forward_ch = Channel.fromPath(params.forward_matches, checkIfExists: true)
    query_ch = Channel.fromPath(params.query_proteins, checkIfExists: true)

    out = FULL_RECIPROCAL_DIAMOND(nr_db_ch, rvdb_db_ch, loci_ch, forward_ch, query_ch)

    out.reciprocal_nr_matches_ch
        .map { it.getName() }
        .collectFile(name: 'reciprocal_nr_manifest.txt', newLine: true, storeDir: params.outdir)

    out.reciprocal_rvdb_matches_ch
        .map { it.getName() }
        .collectFile(name: 'reciprocal_rvdb_manifest.txt', newLine: true, storeDir: params.outdir)

    out.mixed_hits
        .map { it.getName() }
        .collectFile(name: 'mixed_hits_manifest.txt', newLine: true, storeDir: params.outdir)
}
