// Workflow specific parameters

params.diamond_forks = "4"

// Check for outdir & input files

if (file("$params.outdir").exists()) {
    error("Folder '$params.outdir' already exists. Remove or rename it before rerunning, or use --outdir to direct workflow outputs to an alternative directory.")
}

if (!file("$params.assembly_file").exists()) {
    error("Input file '$params.assembly_file' was not found. Either add it to the data directory, or specify another assembly list using --ftp_file")
}

if (!file("$params.reciprocal_nr_db").exists()) {
    error("Input file '$params.reciprocal_nr_db' was not found. Either add it to the data directory, or specify another reciprocal nr database using --reciprocal_nr_db")
}

if (!file("$params.reciprocal_rvdb_db").exists()) {
    error("Input file '$params.reciprocal_rvdb_db' was not found. Either add it to the data directory, or specify another reciprocal RVDB database using --reciprocal_rvdb_db")
}

if (!file("$params.clustered_query_aa").exists()) {
    error("Input file '$params.clustered_query_aa' was not found. Either add it to the data directory, or specify another clustered query file using --clustered_query_aa")
}

if (!file("$params.query_hmms").exists()) {
    error("Input file '$params.query_hmms' was not found. Either add it to the data directory, or specify another query hmm profile using --query_hmms")
}

// Import modules

//include { parse_ftp } from '../modules/parse_ftp.nf'
//include { download_assemblies } from '../modules/download_assemblies.nf'
//include { get_assembly_metadata } from '../modules/get_assembly_metadata.nf'
include { assembly_stats } from '../modules/assembly_stats.nf'
include { diamond } from '../modules/diamond.nf'
include { intersect_domains_merge_extract } from '../modules/intersect_domains_merge_extract.nf'
include { orf_extract } from '../modules/orf_extract.nf'
include { reciprocal_diamond } from '../modules/reciprocal_diamond.nf'
include { genewise } from '../modules/genewise.nf'
//include { build_taxonomy_table } from '../modules/build_taxonomy_table.nf'
include { build_taxonomy_table_reciprocal_only } from '../modules/build_taxonomy_table_reciprocal_only.nf'
include { publish } from '../modules/publish.nf'

// Run local workflow

workflow SINGLEFAST {

    // Define channels
    def clustered_query_ch = Channel.fromPath(params.clustered_query_aa)
    def query_hmm_profiles_ch = Channel.fromPath(params.query_hmms)
    def assembly_ch = Channel.fromPath(params.assembly_file)
    def reciprocal_nr_db_ch = Channel.fromPath(params.reciprocal_nr_db)
    def reciprocal_rvdb_db_ch = Channel.fromPath(params.reciprocal_rvdb_db)

    // Unpack user supplied ftp list and begin downloading assemblies
    //fetched_assembly_files = parse_ftp(ftp_ch) | flatten | download_assemblies

    // Get assembly metadata
    //get_assembly_metadata()

    // Get stats about downloaded assembly files
    assembly_stats(assembly_ch).collectFile(name: 'assembly_stats.tsv', newLine: false, storeDir: "${params.outdir}/sql")

    // Run main EVE search, annotate potential domains, and extract FASTAs
    diamond_out = diamond(assembly_ch.combine(clustered_query_ch))
    intersect_domains_merge_extract(diamond_out.combine(query_hmm_profiles_ch))
    strict_fastas_collected = intersect_domains_merge_extract.out.strict_fa_ch.collect()
    context_fastas_collected = intersect_domains_merge_extract.out.context_fa_ch.collect()
    locus_assembly_map_collected = intersect_domains_merge_extract.out.locus_assembly_map_ch.collect()
    annotated_hits_collected = intersect_domains_merge_extract.out.annot_tsv_ch.collect()

    // Extract ORFs that overlap DIAMOND hits (extending into flanks)
    orfs_collected = orf_extract(intersect_domains_merge_extract.out.context_fa_ch, \
                intersect_domains_merge_extract.out.strict_coords_ch).collect()

    // Reciprocal DIAMOND & and prepare for genewise
    reciprocal_diamond(strict_fastas_collected,
                    context_fastas_collected,
                    reciprocal_nr_db_ch,
                    reciprocal_rvdb_db_ch,
                    annotated_hits_collected,
                    clustered_query_ch)

    pair_subsets = reciprocal_diamond.out.pairs_ch.splitText(by: params.pairs_per_task, file: true)
    best_hit_proteins_val = reciprocal_diamond.out.best_hits_fa_ch.collect()
    strict_fastas_val = reciprocal_diamond.out.merged_fa_ch.collect()
    context_fastas_val = reciprocal_diamond.out.context_fa_ch.collect()
    context_coords_val = reciprocal_diamond.out.context_coords_ch.collect()

    // Frameshift and STOP aware reconstruction of EVE sequences
    genewise(pair_subsets,
            best_hit_proteins_val,
            strict_fastas_val,
            context_fastas_val,
            context_coords_val).collectFile(name: 'genewise.tsv', newLine: false, storeDir: "${params.outdir}/sql")

    build_taxonomy_table_reciprocal_only(reciprocal_diamond.out.reciprocal_nr_matches_ch, reciprocal_diamond.out.reciprocal_rvdb_matches_ch)

    // Produce final outputs
    publish(locus_assembly_map_collected, \
            orfs_collected)

}
