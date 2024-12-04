// Workflow specific parameters

// // Check for outdir & input files

// if (file("$params.outdir").exists()) {
//     error("Folder '$params.outdir' already exists. Remove or rename it before rerunning, or use --outdir to direct workflow outputs to an alternative directory.")
// }

// if (!file("$params.query_file_aa").exists()) {
//     error("Input file '$params.query_file_aa' was not found. Either add it to the data directory, or specify another query file using --query_file_aa")
// }

// if (!file("$params.ftp_file").exists()) {
//     error("Input file '$params.ftp_file' was not found. Either add it to the data directory, or specify another assembly list using --ftp_file")
// }

// if (!file("$params.reciprocal_nr_db").exists()) {
//     error("Input file '$params.reciprocal_nr_db' was not found. Either add it to the data directory, or specify another reciprocal nr database using --reciprocal_nr_db")
// }

// if (!file("$params.reciprocal_rvdb_db").exists()) {
//     error("Input file '$params.reciprocal_rvdb_db' was not found. Either add it to the data directory, or specify another reciprocal RVDB database using --reciprocal_rvdb_db")
// }

// if (!file("$params.phmms").exists()) {
//     error("Input file '$params.phmms' was not found. Either add it to the data directory, or specify another phmm database using --phmms")
// }

// if (!params.max_target_seqs) {
//     param.max_target_seqs = 1000
// }

// Import modules

//include { build_db } from '../modules/build_db.nf'
include { cluster_seqs } from '../modules/cluster_seqs.nf'
include { hmmer } from '../modules/hmmer.nf'
include { parse_ftp } from '../modules/parse_ftp.nf'
include { download_assemblies } from '../modules/download_assemblies.nf'
include { get_assembly_metadata } from '../modules/get_assembly_metadata.nf'
include { assembly_stats } from '../modules/assembly_stats.nf'
include { build_diamond_db as build_query} from '../modules/diamond.nf'
include { build_diamond_db as build_reciprocal} from '../modules/diamond.nf'
include { diamond} from '../modules/diamond.nf'
include { extract_seqs_annotate_matches } from '../modules/intersect_domains_merge_extract.nf'
include { orf_extract } from '../modules/orf_extract.nf'
include { reciprocal_diamond } from '../modules/reciprocal_diamond.nf'
include { genewise } from '../modules/genewise.nf'
include { build_taxonomy_table } from '../modules/build_taxonomy_table.nf'
include { publish } from '../modules/publish.nf'

// Run local workflow

workflow {

    // Define channels
    def query_ch = Channel.fromPath(params.query_file_aa, checkIfExists: true)
    def ftp_ch = Channel.fromPath(params.ftp_file, checkIfExists: true)
    // def reciprocal_nr_db_ch = Channel.fromPath(params.reciprocal_nr_db)
    // def reciprocal_rvdb_db_ch = Channel.fromPath(params.reciprocal_rvdb_db)
    println(query_ch.view())
    println(ftp_ch.view())
    
    // If params.cluster_query, cluster sequences
    query_representatives = (params.cluster_query ? cluster_seqs(query_ch) : query_ch) 
    println(query_representatives.view())
    
    // Build DIAMOND database
    vir_db_ch = (params.build_query_db ? 
                    build_query("query", query_representatives): 
                        Channel.fromPath(params.query_diamond_db, checkIfExists: true))
    
    println(vir_db_ch.view())
    
    // Unpack user supplied ftp list and begin downloading assemblies
    fetched_assembly_files = parse_ftp(ftp_ch) | flatten | download_assemblies
    // println(fetched_assembly_files.view())
    
    // // Get assembly metadata
    // if params.get_assembly_metadata:
    //     get_assembly_metadata()

    // Get stats about downloaded assembly files
    assembly_stats(fetched_assembly_files).collectFile(name: 'assembly_stats.tsv', 
                                                        newLine: false, 
                                                        storeDir: "${params.outdir}/sql")
    
    // Run a DIAMOND using chunks of the genome as query against viral sequences as the database
    diamond_out = diamond(fetched_assembly_files.combine(vir_db_ch))
    println(diamond_out.view())
    //diamond_out = diamond(fetched_assembly_files.combine(build_db.out.vir_db_ch))

    // Extract nucleotide locus, and flanking regions
    //intersect_domains_merge_extract(diamond_out.combine(hmmer.out.query_domains_ch))
    extract_seqs_annotate_matches(diamond_out)
    strict_fastas_collected = extract_seqs_annotate_matches.out.strict_fa_ch.collect()
    context_fastas_collected = extract_seqs_annotate_matches.out.context_fa_ch.collect()
    locus_assembly_map_collected = extract_seqs_annotate_matches.out.locus_assembly_map_ch.collect()

    // // TO DO: Add a step to annotate output with hmmer if user specied hmmer
    // // HMMER run on queries
    // if params.phmms:
    //     def profiles_ch = Channel.fromPath(params.phmms, type: 'dir')
    //     hmmer(profiles_ch, query_representatives)
    //     annotated_hits_collected = extract_seqs_annotate_matches.out.annot_tsv_ch.collect()

    // // Extract ORFs that overlap DIAMOND hits (extending into flanks)
    // orfs_collected = orf_extract(extract_seqs_annotate_matches.out.context_fa_ch, \
    //             extract_seqs_annotate_matches.out.strict_coords_ch).collect()

    // // Reciprocal DIAMOND & and prepare for genewise
    if (params.reciprocal) {

        reciprocal_db = (params.build_reciprocal ? 
                            build_reciprocal("reciprocal", Channel.fromPath(params.reciprocal_db)):
                            Channel.fromPath(params.reciprocal_db, checkIfExists: true))
        
        reciprocal_diamond(strict_fastas_collected,
                        context_fastas_collected,
                        reciprocal_db)

    }
    // if params.reciprocal_diamond:
    //     // For reciprocal databases for databases provided:
    //     reciprocal_db = build_diamond_db(params.reciprocal_seqs).out.vir_db_ch

    // pair_subsets = reciprocal_diamond.out.pairs_ch.splitText(by: params.pairs_per_task, file: true)
    // best_hit_proteins_val = reciprocal_diamond.out.best_hits_fa_ch.collect()
    // strict_fastas_val = reciprocal_diamond.out.merged_fa_ch.collect()
    // context_fastas_val = reciprocal_diamond.out.context_fa_ch.collect()
    // context_coords_val = reciprocal_diamond.out.context_coords_ch.collect()

    // // Frameshift and STOP aware reconstruction of EVE sequences
    // genewise(pair_subsets,
    //         best_hit_proteins_val,
    //         strict_fastas_val,
    //         context_fastas_val,
    //         context_coords_val).collectFile(name: 'genewise.tsv', 
    //                                         newLine: false, 
    //                                         storeDir: "${params.outdir}/sql")

    // // Produce taxonomy table for reciprocal searches and host assemblies
    // build_taxonomy_table(ftp_ch,
    //         get_assembly_metadata.out.assembly_metadata_ch,
    //         reciprocal_diamond.out.reciprocal_nr_matches_ch,
    //         reciprocal_diamond.out.reciprocal_rvdb_matches_ch)

    // // Produce final outputs
    // publish(locus_assembly_map_collected, \
    //         orfs_collected)

}
