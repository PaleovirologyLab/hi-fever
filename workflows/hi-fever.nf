// Workflow specific parameters

params.diamond_forks = "4"

// Check for outdir & input files

if (file("$params.outdir").exists()) {
    error("Folder '$params.outdir' already exists. Remove or rename it before rerunning, or use --outdir to direct workflow outputs to an alternative directory.")
}

if (!file("$params.query_file_aa").exists()) {
    error("Input file '$params.query_file_aa' was not found. Either add it, or specify another query file using --query_file_aa")
}

if (!file("$params.ftp_file").exists()) {
    error("Input file '$params.ftp_file' was not found. Either add it, or specify another assembly list using --ftp_file")
}

if (!file("$params.reciprocal_db").exists()) {
    error("Input file '$params.reciprocal_db' was not found. Either add it, or specify another reciprocal database using --reciprocal_db")
}

if (!file("$params.phmms").exists()) {
    error("Input file '$params.phmms' was not found. Either add it, or specify another phmm database using --phmms")
}

// Import modules

include { build_db } from '../modules/build_db.nf'
include { hmmer } from '../modules/hmmer.nf'
include { parse_ftp } from '../modules/parse_ftp.nf'
include { download_assemblies } from '../modules/download_assemblies.nf'
include { assembly_stats } from '../modules/assembly_stats.nf'
include { diamond } from '../modules/diamond.nf'
include { intersect_domains_merge_extract } from '../modules/intersect_domains_merge_extract.nf'
include { orf_extract } from '../modules/orf_extract.nf'
include { genewise } from '../modules/genewise.nf'
include { reciprocal_diamond } from '../modules/reciprocal_diamond.nf'
include { attempt_genewise_improvement } from '../modules/attempt_genewise_improvement.nf'

// Run local workflow

workflow HIFEVER {

    // Define channels
    def query_ch = Channel.fromPath(params.query_file_aa)
    def profiles_ch = Channel.fromPath(params.phmms, type: 'dir')
    def ftp_ch = Channel.fromPath(params.ftp_file)
    def reciprocal_db_ch = Channel.fromPath(params.reciprocal_db)
    clustered_proteins = Channel.value("$PWD/${params.outdir}/virusdb/DB_clu_rep.fasta")

    // Build clustered DIAMOND query database from user supplied protein FASTA
    build_db(query_ch)

    // HMMER run on clustered queries
    hmmer (profiles_ch, build_db.out.clust_ch)

    // Unpack user supplied ftp list and begin downloading assemblies
    fetched_assembly_files = parse_ftp(ftp_ch) | flatten | download_assemblies

    // Get stats about downloaded assembly files
    assembly_stats(fetched_assembly_files).collectFile(name: 'assembly_stats.txt', newLine: false, storeDir: "$params.outdir")

    // Run main EVE search, annotate potential domains, and extract FASTAs
    diamond_out = diamond(fetched_assembly_files.combine(build_db.out.vir_db_ch))
    intersect_domains_merge_extract(diamond_out.combine(hmmer.out.query_domains_ch))
    strict_fastas_collected = intersect_domains_merge_extract.out.strict_fa_ch.collect()
    context_fastas_collected = intersect_domains_merge_extract.out.context_fa_ch.collect()

    // Extract ORFs that overlap DIAMOND hits (extending into flanks)
    orf_extract (intersect_domains_merge_extract.out.context_fa_ch, \
                intersect_domains_merge_extract.out.strict_coords_ch)

    // Frameshift and STOP aware reconstruction of EVE sequences
    collected_genewise = genewise (intersect_domains_merge_extract.out.annot_tsv_ch, \
                                intersect_domains_merge_extract.out.strict_fa_ch, \
                                intersect_domains_merge_extract.out.context_fa_ch, \
                                intersect_domains_merge_extract.out.context_coords_ch, \
                                clustered_proteins) | \
                                collect

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
