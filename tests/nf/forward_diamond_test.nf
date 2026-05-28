nextflow.enable.dsl=2

include { FORWARD_DIAMOND } from '../../modules/forward_diamond.nf'

workflow {
    assembly_ch = Channel.fromPath(params.assembly, checkIfExists: true)
    db_ch = Channel.fromPath(params.query_db, checkIfExists: true)

    meta_ch = assembly_ch.map { asm ->
        def name = asm.name
        name = name.replaceFirst(/\.gz$/, '')
        name = name.replaceFirst(/\.(fa|fna|fasta)$/, '')
        def meta = [id: name]
        return [meta, asm]
    }

    out = FORWARD_DIAMOND(meta_ch.combine(db_ch))

    out
        .map { meta, file -> file.getName() }
        .collectFile(name: 'forward_manifest.txt', newLine: true, storeDir: params.outdir)
}
