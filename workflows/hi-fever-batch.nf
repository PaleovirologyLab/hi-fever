// Workflow specific parameters

params.diamond_forks = "5"
params.bucket_name = "hifeverbucket"

// Import modules

include { build_db } from '../modules/xgc_build_db.nf'
include { hmmer } from '../modules/xgc_hmmer.nf'
include { parse_ftp } from '../modules/xgc_parse_ftp.nf'
include { download_assemblies } from '../modules/xgc_download_assemblies.nf'
include { assembly_stats } from '../modules/xgc_assembly_stats.nf'
include { diamond } from '../modules/xgc_diamond.nf'
include { intersect_domains_merge_extract } from '../modules/xgc_intersect_domains_merge_extract.nf'
include { orf_extract } from '../modules/xgc_orf_extract.nf'
include { genewise } from '../modules/xgc_genewise.nf'
include { reciprocal_diamond } from '../modules/xgc_reciprocal_diamond.nf'
include { attempt_genewise_improvement } from '../modules/xgc_attempt_genewise_improvement.nf'

// Run cloud batch workflow

workflow BATCH {

    // Check that outdir does not exist. If not, proceed. If yes, warn.
    if (!file("$params.outdir").exists()) {

        // Define channels
        def query_ch = Channel.fromPath(params.query_file_aa)
        def profiles_ch = Channel.fromPath(params.phmms, type: 'dir')
        def ftp_ch = Channel.fromPath(params.ftp_file)
        def reciprocal_db_ch = Channel.fromPath(params.reciprocal_db)
        clustered_proteins = Channel.value("gs://${params.bucket_name}/${params.outdir}/virusdb/DB_clu_rep.fasta")

        // Build clustered DIAMOND query database from user supplied protein FASTA
        build_db(query_ch)

        // HMMER run on clustered queries
        hmmer (profiles_ch, build_db.out.clust_ch)

        // Unpack user supplied ftp list and begin downloading assemblies
        fetched_assembly_files = parse_ftp(ftp_ch) | flatten | download_assemblies

        // Get stats about downloaded assembly files
        assembly_stats(fetched_assembly_files).collectFile(name: 'assembly_stats.txt', newLine: false, storeDir: "gs://${params.bucket_name}/${params.outdir}")

        // Run main EVE search, annotate potential domains, and extract FASTAs
        diamond_out = diamond(fetched_assembly_files.combine(build_db.out.vir_db_ch))
        intersect_domains_merge_extract(diamond_out.combine(hmmer.out.query_domains_ch))
        strict_fastas_collected = intersect_domains_merge_extract.out.strict_fa_ch.collect()
        context_fastas_collected = intersect_domains_merge_extract.out.context_fa_ch.collect()

        // Extract ORFs that overlap DIAMOND hits (extending into flanks)
        orf_extract (intersect_domains_merge_extract.out.context_fa_ch, \
                     intersect_domains_merge_extract.out.strict_coords_ch)

        // Frameshift and STOP aware reconstruction of EVE sequences

        collected_genewise = genewise (intersect_domains_merge_extract.out.idme_tuple.combine(clustered_proteins)) | collect

        // Reciprocal DIAMOND
        reciprocal_diamond (collected_genewise, reciprocal_db_ch)

        // Attempt genewise improvement if reciprocal best hit changed
        attempt_genewise_improvement (reciprocal_diamond.out.original_genewise_ch, \
                                      reciprocal_diamond.out.reciprocal_hits_ch, \
                                      reciprocal_diamond.out.reciprocal_fasta_ch, \
                                      reciprocal_db_ch, \
                                      strict_fastas_collected, \
                                      context_fastas_collected)

    }

    else {

        println ("Folder '$params.outdir' already exists. Remove or rename it before rerunning, or use --outdir to direct workflow outputs to an alternative directory.")

    }

}
