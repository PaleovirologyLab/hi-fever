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
include { get_assembly_metadata; get_assembly_metadata_from_accn} from '../modules/get_assembly_metadata.nf'
include { assembly_stats } from '../modules/assembly_stats.nf'
// DIAMOND processes
include { build_diamond_db as build_query} from '../modules/diamond.nf'
include { build_diamond_db as build_reciprocal} from '../modules/diamond.nf'
include { diamond } from '../modules/diamond.nf'
include { reciprocal_diamond_single } from '../modules/diamond.nf'
include { reciprocal_diamond_double } from '../modules/diamond.nf'
// Concatenating and extracting best hits when reciprocal DIAMOND
include { merge_seqs_loci} from '../modules/post_processing_files.nf'
include {find_best_hits} from '../modules/post_processing_files.nf'
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
    // // def reciprocal_nr_db_ch = Channel.fromPath(params.reciprocal_nr_db)
    // // def reciprocal_rvdb_db_ch = Channel.fromPath(params.reciprocal_rvdb_db)
    // println(query_ch.view())
    // println(ftp_ch.view())
    
    // // If params.cluster_query, cluster sequences
    query_representatives = (params.cluster_query ? cluster_seqs(query_ch) : query_ch) 
    // println(query_representatives.view())
    
    // Build DIAMOND database
    vir_db_ch = (params.build_query_db ? 
                    build_query("query", query_representatives): 
                        Channel.fromPath(params.query_diamond_db, checkIfExists: true))
    
    // println(vir_db_ch.view())
    
    // // Unpack user supplied ftp list and begin downloading assemblies
    fetched_assembly_files = parse_ftp(ftp_ch) | flatten | download_assemblies
    // fetched_assembly_files = Channel.fromPath("/home/biology/biol0273/Projects/hi-fever/test_data/GCF_035594765.1_rAnoCar3.1.pri_genomic.fna.gz")
    println(fetched_assembly_files.view())
    
    // // // Get assembly metadata
    // // if params.get_assembly_metadata:
    // //     get_assembly_metadata()

    // // Get stats about downloaded assembly files
    assembly_stats = assembly_stats(fetched_assembly_files).collectFile(name: 'assembly_stats.tsv', 
                                                        newLine: false, 
                                                        storeDir: "${params.outdir}/sql")

    get_assembly_metadata_from_accn(assembly_stats).collectFile(name: 'assembly_metadata.tsv', 
                                                        newLine: false, 
                                                        storeDir: "${params.outdir}/sql")

    
    // Run a DIAMOND using chunks of the genome as query against viral sequences as the database
    diamond_out = diamond(fetched_assembly_files.combine(vir_db_ch))
    
    // Extract nucleotide locus, and flanking regions
    extract_seqs_annotate_matches(diamond_out)

    // Inputs for reciprocal DIAMOND:
    forward_matches_collected = extract_seqs_annotate_matches.out.forward_matches.collect()
    strict_fastas_collected = extract_seqs_annotate_matches.out.strict_fa_ch.collect()
    context_fastas_collected = extract_seqs_annotate_matches.out.context_fa_ch.collect()
    
    // Files to publish:
    locus_assembly_map_collected = extract_seqs_annotate_matches.out.locus_assembly_map_ch.collect()
    strict_coords_ch = extract_seqs_annotate_matches.out.strict_coords_ch.collect()

    // forward_matches_collected = Channel.fromPath("/home/biology/biol0273/Projects/hi-fever/test_data/extracted/*_forward_matches.dmnd.annot.tsv")
    // strict_fastas_collected = Channel.fromPath("/home/biology/biol0273/Projects/hi-fever/test_data/extracted/*_strict.fasta")
    // context_fastas_collected = Channel.fromPath("/home/biology/biol0273/Projects/hi-fever/test_data/extracted/*_context.fasta")
    // locus_assembly_map_collected = Channel.fromPath("/home/biology/biol0273/Projects/hi-fever/test_data/extracted/*_locus_assembly_map.tsv")
    // strict_coords_ch = Channel.fromPath("/home/biology/biol0273/Projects/hi-fever/test_data/extracted/strict_coords.bed")

    // TO DO: Add a step to annotate output with hmmer if user specied hmmer
    // // HMMER run on queries
    // if params.phmms:
    //     def profiles_ch = Channel.fromPath(params.phmms, type: 'dir')
    //     hmmer(profiles_ch, query_representatives)
    //     forward_matches_collected = extract_seqs_annotate_matches.out.annot_tsv_ch.collect()

    // // Extract ORFs that overlap DIAMOND hits (extending into flanks)
    // orfs_collected = orf_extract(extract_seqs_annotate_matches.out.context_fa_ch, \
    //             extract_seqs_annotate_matches.out.strict_coords_ch).collect()

    merge_seqs_loci(strict_fastas_collected, context_fastas_collected)
    loci_merged_fa = merge_seqs_loci.out.loci_merged_fa.collect()
    loci_merged_context_gz = merge_seqs_loci.out.loci_merged_context_gz.collect()
    all_context_coords_bed = merge_seqs_loci.out.all_context_coords_bed.collect()

    // Reciprocal DIAMOND & and prepare for genewise
    if (params.reciprocal) {

        reciprocal_db = (params.build_reciprocal ? 
                            build_reciprocal("reciprocal", Channel.fromPath(params.reciprocal_db)):
                            Channel.fromPath(params.reciprocal_db, checkIfExists: true))
        println(reciprocal_db.view())
        
        reciprocal_diamond_single(reciprocal_db, loci_merged_fa)

        reciprocal_matches = reciprocal_diamond_single.out.reciprocal_matches.collect()
        reciprocal_seqs = reciprocal_diamond_single.out.reciprocal_seqs.collect()
        reciprocal_hits = reciprocal_diamond_single.out.reciprocal_hits.collect()
        // loci_merged_fa_gz = reciprocal_diamond_single.out.loci_merged_fa_gz.collect()
        find_best_hits(forward_matches_collected, 
                        query_representatives, reciprocal_matches, 
                        reciprocal_seqs, reciprocal_hits)
    }


    pair_subsets = find_best_hits.out.best_pairs_txt.splitText(by: params.pairs_per_task, file: true)
    best_hit_proteins_val = find_best_hits.out.best_hits_fa_ch.collect()
    
    // Frameshift and STOP aware reconstruction of EVE sequences
    genewise(pair_subsets,
            best_hit_proteins_val,
            loci_merged_fa,
            loci_merged_context_gz,
            all_context_coords_bed).collectFile(name: 'genewise.tsv', 
                                            newLine: false, 
                                            storeDir: "${params.outdir}/sql")

    // Produce taxonomy table for reciprocal searches and host assemblies
    // build_taxonomy_table(ftp_ch,
    //         get_assembly_metadata.out.assembly_metadata_ch,
    //         reciprocal_diamond.out.reciprocal_nr_matches_ch,
    //         reciprocal_diamond.out.reciprocal_rvdb_matches_ch)

    // // Produce final outputs
    // publish(locus_assembly_map_collected, \
    //         orfs_collected)

}
