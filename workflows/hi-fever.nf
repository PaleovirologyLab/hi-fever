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

// Process file for outputs and publishing
include { build_taxonomy_table } from '../modules/build_taxonomy_table.nf'
include { build_diamond_taxonomy_table } from '../modules/build_taxonomy_table.nf'
include { build_host_lineage_table } from '../modules/build_taxonomy_table.nf'

// Publish concatenated tables with results of all hosts
include { concatenate_publish_tables as pulish_predicted_orfs} from '../modules/post_processing_files.nf'
include { concatenate_publish_tables as pulish_forward_diamond} from '../modules/post_processing_files.nf'
include { concatenate_publish_tables as pulish_assembly_map} from '../modules/post_processing_files.nf'

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
        metadata_channel = get_metadata(assembly_stats)
        build_host_lineage_table(metadata_channel.assembly_metadata)        

    }  
    else if (params.get_all_metada) {
        // Download assembly metadata for all eukaryots, subset target
        get_assembly_metadata_all()

    }
    
    // Run a DIAMOND using chunks of the genome as query against and viral sequences as database
    forward_diamond_out = forward_diamond(fetched_assembly_files.combine(vir_db_ch))
    
    // Extract genome locus with hits, and their flanking regions
    extract_seqs_outputs = extract_seqs_annotate_matches(forward_diamond_out)

    // Inputs for reciprocal DIAMOND:
    forward_matches = extract_seqs_outputs.forward_matches.collect()
    strict_fastas_collected = extract_seqs_outputs.strict_fa_ch.collect()
    context_fastas_collected = extract_seqs_outputs.context_fa_ch.collect()
    
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
        
        // run reciprocal DIAMOND and publish results
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
        all_diamond_hits = find_best_diamond_hits.out.forward_plus_reciprocal_dmd_hits.collect()
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
        all_diamond_hits = reciprocal_diamond_full.out.mixed_hits.collect()
        
    }
    
    // Reconstrunction of encoded proteins. Returns STOP codons, frameshifts and indels
    genewise(best_pairs_subsets,
            best_hit_proteins_val,
            loci_merged_fa,
            loci_merged_context_gz,
            all_context_coords_bed).collectFile(name: 'genewise.tsv', 
                                            newLine: false, 
                                            storeDir: "${params.outdir}/sql")

    // Make taxonomy and publish table for proteins
    build_diamond_taxonomy_table(all_diamond_hits)

    // Publish files
    cat_forward = pulish_forward_diamond(forward_matches, "forward-matches.dmnd.annot.tsv")
    predicted_orfs = orf_extract(extract_seqs_outputs.context_fa_ch, 
                                 extract_seqs_outputs.strict_coords_ch)
    cat_orfs = pulish_predicted_orfs(predicted_orfs.orfs.collect(), "predicted_orfs.tsv")
    locus_assembly_maps = extract_seqs_outputs.locus_assembly_map_ch.collect()
    cat_assembly_map = pulish_assembly_map(locus_assembly_maps, "locus_assembly_map.tsv")

}
