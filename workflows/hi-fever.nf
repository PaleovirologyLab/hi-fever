// Parse inputs
include { parse_ftp } from '../modules/process_host_info.nf'
include { check_file_type as check_query_type } from '../modules/utils.nf'
include { check_file_type as check_reciprocal_type} from '../modules/utils.nf'
include { cluster_seqs } from '../modules/cluster_seqs.nf'

// Process host genomes and get their taxonomical information
include { download_assemblies } from '../modules/process_host_info.nf'
include { download_extract_host_metadata; get_metadata} from '../modules/process_host_info.nf'
include { assembly_stats } from '../modules/process_host_info.nf'
include { fetch_host_taxonomy } from '../modules/process_host_info.nf'
include { build_host_taxonomy_table } from '../modules/process_host_info.nf'

// DIAMOND-related process
include { build_diamond_db as build_query} from '../modules/diamond.nf'
include { build_diamond_db as build_reciprocal} from '../modules/diamond.nf'
include { forward_diamond } from '../modules/diamond.nf'
include { single_reciprocal_diamond; full_reciprocal_diamond } from '../modules/diamond.nf'
include { find_best_diamond_hits } from '../modules/diamond.nf'

// Concatenate and extract best hits after forward DIAMOND
include { extract_seqs_annotate_matches } from '../modules/intersect_domains_merge_extract.nf'
include { merge_seqs_loci} from '../modules/merge_seqs_loci.nf'

// Build ancestral protein
include { genewise } from '../modules/genewise.nf'

// Extract orfs from forward DIAMOND
include { orf_extract } from '../modules/orf_extract.nf'

// Getting hits taxonomy tables
include { build_hits_taxonomy_table } from '../modules/hits_taxonomy.nf'
include { fetch_hits_taxonomy_from_accns } from '../modules/hits_taxonomy.nf'

// Publish concatenated tables with results of all hosts
include { concatenate_publish_tables as publish_predicted_orfs} from '../modules/utils.nf'
include { concatenate_publish_tables as publish_forward_diamond} from '../modules/utils.nf'
include { concatenate_publish_tables as publish_assembly_map} from '../modules/utils.nf'

// Run local workflow

workflow HIFEVER {

    // Define channels
    def query_ch = Channel.fromPath(params.query_file_aa, checkIfExists: true)
    def ftp_ch = Channel.fromPath(params.ftp_file, checkIfExists: true)

    // If params.cluster_query, cluster sequences
    query_proteins = (params.cluster_query ? cluster_seqs(query_ch) : query_ch) 

    // Build DIAMOND database with queries if fasta file
    query_type = query_ch | check_query_type
    vir_db_ch = (query_type == 'fasta' ? 
                    build_query("query", query_proteins): query_proteins)
    
    // Unpack ftp list, download assemblies
    fetched_assembly_files = parse_ftp(ftp_ch) | flatten | download_assemblies
   
    // Get stats about downloaded assembly files
    assembly_stats = assembly_stats(fetched_assembly_files).collectFile(name: 'assembly_stats.tsv', 
                                                        newLine: false, 
                                                        storeDir: "${params.outdir}/sql")

    if (!params.get_all_metadata) {
        // Download metadata only for genomes on ftp file
        metadata_channel = get_metadata(assembly_stats)
        fetch_host_taxonomy(metadata_channel.assembly_metadata)        

    }  
    else {
        // Download assembly metadata for all eukaryots
        download_extract_host_metadata()
        def ncbi_tax_table = Channel.fromPath(params.ncbi_taxonomy_table, checkIfExists: true)
        build_host_taxonomy_table( ftp_ch, 
                                   download_extract_host_metadata.out.assembly_metadata_ch, 
                                   ncbi_tax_table)

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
    if (params.custom_reciprocal) {
        
        def reciprocal_ch = Channel.fromPath(params.custom_reciprocal_db, checkIfExists: true)
        reciprocal_type = reciprocal_ch | check_reciprocal_type

        // Create reciprocal database if input is in fasta, else use reciprocal_db
        reciprocal_db = (reciprocal_type == 'fasta' ? 
                            build_reciprocal("reciprocal", reciprocal_ch): reciprocal_ch)
        
        // run reciprocal DIAMOND and publish results
        single_reciprocal_diamond(reciprocal_db, loci_merged_fa)    
        reciprocal_matches = single_reciprocal_diamond.out.reciprocal_matches.collect()
        reciprocal_seqs = single_reciprocal_diamond.out.reciprocal_seqs.collect()
        reciprocal_hits = single_reciprocal_diamond.out.reciprocal_hits.collect()

        // Find best hits    
        find_best_diamond_hits(forward_matches, 
                        query_proteins, reciprocal_matches, 
                        reciprocal_seqs, reciprocal_hits)

        best_pairs_subsets = find_best_diamond_hits.out.best_pairs_txt.splitText(
                                                            by: params.pairs_per_task, 
                                                            file: true)

        best_hit_proteins_val = find_best_diamond_hits.out.best_hits_fa_ch.collect()
        all_diamond_hits = find_best_diamond_hits.out.forward_plus_reciprocal_dmnd_hits.collect()
        
        // Make taxonomy and publish table for proteins
        hits_taxonomy = fetch_hits_taxonomy_from_accns(all_diamond_hits)
    }
    else {
        // Reciprocal DIAMOND with rvdb and nr protein databases
        def reciprocal_nr_db_ch = Channel.fromPath(params.reciprocal_nr_db, checkIfExists: true)
        def reciprocal_rvdb_db_ch = Channel.fromPath(params.reciprocal_rvdb_db, checkIfExists: true)
        
        full_reciprocal_diamond(reciprocal_nr_db_ch, reciprocal_rvdb_db_ch,
                                loci_merged_fa, forward_matches,
                                query_proteins)
        best_pairs_subsets = full_reciprocal_diamond.out.best_pairs_txt.splitText(
                                                            by: params.pairs_per_task, 
                                                            file: true)
        best_hit_proteins_val = full_reciprocal_diamond.out.best_hits_fa_ch.collect()
        all_diamond_hits = full_reciprocal_diamond.out.mixed_hits.collect()

        // Read taxonomy table to build hits taxonomy
        def ncbi_tax_table_hits = Channel.fromPath(params.ncbi_taxonomy_table, checkIfExists: true)
        
        // Build hits taxonomy from annotated diamond database
        hits_taxonomy = build_hits_taxonomy_table( full_reciprocal_diamond.out.reciprocal_nr_matches_ch, 
                                                   full_reciprocal_diamond.out.reciprocal_rvdb_matches_ch, 
                                                   ncbi_tax_table_hits)
        
    }
    
    // Reconstruction of encoded proteins. Returns number of STOP codons, frameshifts and indels
    genewise(best_pairs_subsets,
            best_hit_proteins_val,
            loci_merged_fa,
            loci_merged_context_gz,
            all_context_coords_bed).collectFile(name: 'genewise.tsv', 
                                            newLine: false, 
                                            storeDir: "${params.outdir}/sql")

    // Publish files
    cat_forward = publish_forward_diamond(forward_matches, "forward-matches.dmnd.annot.tsv")
    predicted_orfs = orf_extract(extract_seqs_outputs.context_fa_ch, 
                                 extract_seqs_outputs.strict_coords_ch)

    cat_orfs = publish_predicted_orfs(predicted_orfs.orfs.collect(), "predicted_orfs.tsv")
    locus_assembly_maps = extract_seqs_outputs.locus_assembly_map_ch.collect()
    cat_assembly_map = publish_assembly_map(locus_assembly_maps, "locus_assembly_map.tsv")

}
