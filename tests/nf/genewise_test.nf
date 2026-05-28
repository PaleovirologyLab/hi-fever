nextflow.enable.dsl=2

include { GENEWISE } from '../../modules/genewise.nf'

workflow {
    pair_ch = Channel.fromPath(params.pair_subsets, checkIfExists: true)
    proteins_ch = Channel.fromPath(params.best_hit_proteins, checkIfExists: true)
    strict_ch = Channel.fromPath(params.strict_fastas, checkIfExists: true)
    context_ch = Channel.fromPath(params.context_fastas, checkIfExists: true)
    coords_ch = Channel.fromPath(params.context_coords, checkIfExists: true)

    out = GENEWISE(pair_ch, proteins_ch, strict_ch, context_ch, coords_ch)

    out.genewise_file
        .map { it.getName() }
        .collectFile(name: 'genewise_manifest.txt', newLine: true, storeDir: params.outdir)
}
