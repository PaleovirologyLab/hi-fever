nextflow.enable.dsl=2

include { EXTRACT_SEQS_ANNOTATE_MATCHES } from '../../modules/intersect_domains_merge_extract.nf'

workflow {
    assembly_db = file(params.assembly_db)
    diamond_tsv = file(params.diamond_tsv)
    meta_id = params.meta_id ?: "assembly"
    meta = [id: meta_id]

    input_ch = Channel.of([meta, diamond_tsv, assembly_db])
    out = EXTRACT_SEQS_ANNOTATE_MATCHES(input_ch)

    out.forward_matches
        .map { it.getName() }
        .collectFile(name: 'forward_annot_manifest.txt', newLine: true, storeDir: params.outdir)

    out.strict_fa_ch
        .map { it.getName() }
        .collectFile(name: 'strict_manifest.txt', newLine: true, storeDir: params.outdir)

    out.context_fa_ch
        .map { it.getName() }
        .collectFile(name: 'context_manifest.txt', newLine: true, storeDir: params.outdir)

    out.locus_assembly_map_ch
        .map { it.getName() }
        .collectFile(name: 'locus_map_manifest.txt', newLine: true, storeDir: params.outdir)

    out.strict_coords_ch
        .map { it.getName() }
        .collectFile(name: 'strict_coords_manifest.txt', newLine: true, storeDir: params.outdir)
}
