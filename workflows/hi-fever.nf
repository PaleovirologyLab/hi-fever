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

// INPUTS processes
include { cluster_seqs } from '../modules/cluster_seqs.nf'
// include { hmmer } from '../modules/hmmer.nf'
include { parse_ftp } from '../modules/parse_ftp.nf'
include { download_assemblies } from '../modules/download_assemblies.nf'
include { get_assembly_metadata_all; get_metadata} from '../modules/get_assembly_metadata.nf'
include { assembly_stats } from '../modules/assembly_stats.nf'

// DIAMOND processes
include { build_diamond_db as build_query} from '../modules/diamond.nf'
include { build_diamond_db as build_reciprocal} from '../modules/diamond.nf'
include { forward_diamond } from '../modules/diamond.nf'
include { reciprocal_diamond; reciprocal_diamond_full } from '../modules/diamond.nf'
include { find_best_diamond_hits } from '../modules/diamond.nf'

// Concatenating and extracting best hits when reciprocal DIAMOND
include { merge_seqs_loci} from '../modules/post_processing_files.nf'
include { extract_seqs_annotate_matches } from '../modules/intersect_domains_merge_extract.nf'
include { orf_extract } from '../modules/orf_extract.nf'
include { genewise } from '../modules/genewise.nf'
include { build_taxonomy_table } from '../modules/build_taxonomy_table.nf'
include { publish } from '../modules/publish.nf'

// Run local workflow

workflow HIFEVER {

    // Define channels
    def query_ch = Channel.fromPath(params.query_file_aa, checkIfExists: true)
    def ftp_ch = Channel.fromPath(params.ftp_file, checkIfExists: true)

    // If params.cluster_query, cluster sequences
    query_proteins = (params.cluster_query ? cluster_seqs(query_ch) : query_ch) 
    
    // Build DIAMOND database
    vir_db_ch = (params.build_query_db ? 
                    build_query("query", query_proteins): 
                        Channel.fromPath(params.query_diamond_db, checkIfExists: true))
    
    
    // Unpack ftp list, download assemblies
    fetched_assembly_files = parse_ftp(ftp_ch) | flatten | download_assemblies
   
    // Get stats about downloaded assembly files
    assembly_stats = assembly_stats(fetched_assembly_files).collectFile(name: 'assembly_stats.tsv', 
                                                        newLine: false, 
                                                        storeDir: "${params.outdir}/sql")

    // Get assembly metadata
    if (!params.dont_get_metada) {
        // For target genomes, requires email
        get_metadata(assembly_stats).collectFile(name: 'assembly_metadata.tsv', 
                                                 newLine: false, 
                                                 storeDir: "${params.outdir}/sql")

    }  
    else if (params.get_all_metada) {
        // Download assembly metadata for all eukaryots, subset target
        get_assembly_metadata_all()

    }
    
    // Run a DIAMOND using chunks of the genome as query against and viral sequences as database
    forward_diamond_out = forward_diamond(fetched_assembly_files.combine(vir_db_ch))
    
    // Extract genome locus with hits, and their flanking regions
    extract_seqs_outputs = extract_seqs_annotate_matches(forward_diamond_out)

    orfs_collected = orf_extract(extract_seqs_outputs.context_fa_ch, 
                                 extract_seqs_outputs.strict_coords_ch)

    // Inputs for reciprocal DIAMOND:
    forward_matches = extract_seqs_outputs.forward_matches.collect()
    strict_fastas_collected = extract_seqs_outputs.strict_fa_ch.collect()
    context_fastas_collected = extract_seqs_outputs.context_fa_ch.collect()
    
    // Files to publish
    locus_assembly_map_collected = extract_seqs_outputs.locus_assembly_map_ch.collect()

    // Concatenate sequences into single file
    merge_seqs_loci(strict_fastas_collected, context_fastas_collected)
    loci_merged_fa = merge_seqs_loci.out.loci_merged_fa.collect()
    loci_merged_context_gz = merge_seqs_loci.out.loci_merged_context_gz.collect()
    all_context_coords_bed = merge_seqs_loci.out.all_context_coords_bed.collect()

    // Run reciprocal DIAMOND
    if (!params.full_reciprocal) {
        
        // build reciprocal if provided by user, else, use default provided by us
        reciprocal_db = (!params.dont_build_reciprocal ? 
                            build_reciprocal("reciprocal", Channel.fromPath(params.reciprocal_db)):
                            Channel.fromPath(params.reciprocal_db, checkIfExists: true))
        
        // run reciprocal DIAMOND
        reciprocal_diamond(reciprocal_db, loci_merged_fa)    
        reciprocal_matches = reciprocal_diamond.out.reciprocal_matches.collect()
        reciprocal_seqs = reciprocal_diamond.out.reciprocal_seqs.collect()
        reciprocal_hits = reciprocal_diamond.out.reciprocal_hits.collect()

        // Find best hits    
        find_best_diamond_hits(forward_matches, 
                        query_proteins, reciprocal_matches, 
                        reciprocal_seqs, reciprocal_hits)

        best_pairs_subsets = find_best_diamond_hits.out.best_pairs_txt.splitText(
                                                            by: params.pairs_per_task, 
                                                            file: true)

        best_hit_proteins_val = find_best_diamond_hits.out.best_hits_fa_ch.collect()
    }
    else {
        // Reciprocal DIAMOND with rvdb and nr protein databases
        def reciprocal_nr_db_ch = Channel.fromPath(params.reciprocal_nr_db)
        def reciprocal_rvdb_db_ch = Channel.fromPath(params.reciprocal_rvdb_db)
        
        reciprocal_diamond_full(reciprocal_nr_db_ch, reciprocal_rvdb_db_ch,
                                loci_merged_fa, forward_matches,
                                query_proteins)
        best_pairs_subsets = reciprocal_diamond_full.out.best_pairs_txt.splitText(
                                                            by: params.pairs_per_task, 
                                                            file: true)
        best_hit_proteins_val = reciprocal_diamond_full.out.best_hits_fa_ch.collect()
        
    }
    
    // Reconstrunction of encoded proteins. Returns STOP codons, frameshifts and indels
    genewise(best_pairs_subsets,
            best_hit_proteins_val,
            loci_merged_fa,
            loci_merged_context_gz,
            all_context_coords_bed).collectFile(name: 'genewise.tsv', 
                                            newLine: false, 
                                            storeDir: "${params.outdir}/sql")

    // TO DO: Add a step to annotate output with hmmer if user specied hmmer
    // // HMMER run on queries
    // if params.phmms:
    //     def profiles_ch = Channel.fromPath(params.phmms, type: 'dir')
    //     hmmer(profiles_ch, query_proteins)
    //     forward_matches = extract_seqs_annotate_matches.out.annot_tsv_ch.collect()


    // Extract ORFs that overlap DIAMOND hits (extending into flanks)
    // orfs_collected = orf_extract(extract_seqs_annotate_matches.out.context_fa_ch, \
    //             extract_seqs_annotate_matches.out.strict_coords_ch).collect()


    // Create predicted_ORFS by looking at which 
    


    // Produce taxonomy table for reciprocal searches and host assemblies
    // build_taxonomy_table(ftp_ch,
    //         get_assembly_metadata.out.assembly_metadata_ch,
    //         reciprocal_diamond.out.reciprocal_nr_matches_ch,
    //         reciprocal_diamond.out.reciprocal_rvdb_matches_ch)

    // // Produce final outputs
    // publish(locus_assembly_map_collected, \
    //         orfs_collected)

}
