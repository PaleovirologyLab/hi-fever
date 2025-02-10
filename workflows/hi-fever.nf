// Parse inputs
include { PARSE_FTP } from '../modules/process_host_info.nf'
include { CLUSTER_SEQS } from '../modules/cluster_seqs.nf'

// Process host genomes and get their taxonomical information
include { DOWNLOAD_ASSEMBLIES } from '../modules/process_host_info.nf'
include { DOWNLOAD_EXTRACT_HOST_METADATA } from '../modules/process_host_info.nf'
include { GET_METADATA } from '../modules/process_host_info.nf'
include { ASSEMBLY_STATS } from '../modules/process_host_info.nf'
include { FETCH_HOST_TAXONOMY } from '../modules/process_host_info.nf'
include { BUILD_HOST_TAXONOMY_TABLE } from '../modules/process_host_info.nf'
include { BLASTDB } from '../modules/process_host_info.nf'

// DIAMOND-related process
include { BUILD_DIAMOND_DB as BUILD_QUERY} from '../modules/diamond.nf'
include { BUILD_DIAMOND_DB as BUILD_RECIPROCAL} from '../modules/diamond.nf'
include { FORWARD_DIAMOND } from '../modules/diamond.nf'
include { SINGLE_RECIPROCAL_DIAMOND } from '../modules/diamond.nf'
include { FULL_RECIPROCAL_DIAMOND } from '../modules/diamond.nf'
include { FIND_BEST_DIAMOND_HITS } from '../modules/diamond.nf'

// Concatenate and extract best hits after forward DIAMOND
include { EXTRACT_SEQS_ANNOTATE_MATCHES } from '../modules/intersect_domains_merge_extract.nf'
include { MERGE_SEQS_LOCI} from '../modules/merge_seqs_loci.nf'

// Build ancestral protein
include { GENEWISE } from '../modules/genewise.nf'

// Extract orfs from forward DIAMOND
include { ORF_EXTRACT } from '../modules/orf_extract.nf'

// Getting hits taxonomy tables
include { BUILD_HITS_TAXONOMY_TABLE } from '../modules/hits_taxonomy.nf'
include { FETCH_HITS_TAXONOMY_FROM_ACCNS } from '../modules/hits_taxonomy.nf'

// Publish concatenated tables with results of all hosts
include { CONCATENATE_PUBLISH_TABLES as PUBLISH_PREDICTED_ORFS} from '../modules/utils.nf'
include { CONCATENATE_PUBLISH_TABLES as PUBLISH_FORWARD_DIAMOND} from '../modules/utils.nf'
include { CONCATENATE_PUBLISH_TABLES as PUBLISH_ASSEMBLY_MAP} from '../modules/utils.nf'

// Run local workflow

workflow HIFEVER {

	// Define channels

		def query_ch = Channel.fromPath("${params.data_path}/${params.query_file_aa}", checkIfExists: true)
		def ftp_ch = Channel.fromPath("${params.data_path}/${params.ftp_file}", checkIfExists: true)

	// If params.cluster_query, cluster sequences and reassign query_ch, else no change

		query_ch = params.cluster_query ? CLUSTER_SEQS(query_ch) : query_ch

	// Build DIAMOND database from queries

		def vir_db_ch = params.query_db ? query_db : BUILD_QUERY("query", query_ch)

    // Unpack ftp list, download assemblies
    fetched_assembly_files = PARSE_FTP(ftp_ch) | flatten | DOWNLOAD_ASSEMBLIES

	// Add assembly metadata
	assembly_with_metadata = fetched_assembly_files.map { assembly ->
									def meta = [
									id: assembly.baseName
									]
									return [meta, assembly]
							}

    // Get stats about downloaded assembly files
    assembly_stats = ASSEMBLY_STATS(fetched_assembly_files).collectFile(name: 'assembly_stats.tsv',
                                                        newLine: false,
                                                        storeDir: "${params.outdir}/sql")

    if (!params.get_all_metadata) {
        // Download metadata only for genomes on ftp file
        metadata_channel = GET_METADATA(assembly_stats)
        FETCH_HOST_TAXONOMY(metadata_channel.assembly_metadata)

    }  
    else {
        // Download assembly metadata for all eukaryots
        DOWNLOAD_EXTRACT_HOST_METADATA()
        def ncbi_tax_table = Channel.fromPath("${params.data_path}/${params.ncbi_taxonomy_table}", checkIfExists: true)
        BUILD_HOST_TAXONOMY_TABLE( ftp_ch,
                                   DOWNLOAD_EXTRACT_HOST_METADATA.out.assembly_metadata_ch,
                                   ncbi_tax_table)

    }

    // Run a DIAMOND using chunks of the genome as query against and viral sequences as database
    forward_diamond_out = FORWARD_DIAMOND(assembly_with_metadata.combine(vir_db_ch))

	// Make BLAST db for assembly
	BLASTDB(assembly_with_metadata)

	// Rejoin dmnd.nsq with nsq db
	rejoinDb = forward_diamond_out.join(BLASTDB.out, by: [0]) // Joins on shared metadata (assembly)


    // Extract genome locus with hits, and their flanking regions
    extract_seqs_outputs = EXTRACT_SEQS_ANNOTATE_MATCHES(rejoinDb)

    // Inputs for reciprocal DIAMOND:
    forward_matches = extract_seqs_outputs.forward_matches.collect()
    strict_fastas_collected = extract_seqs_outputs.strict_fa_ch.collect()
    context_fastas_collected = extract_seqs_outputs.context_fa_ch.collect()
    
    // Concatenate sequences into single file
    MERGE_SEQS_LOCI(strict_fastas_collected, context_fastas_collected)
    loci_merged_fa = MERGE_SEQS_LOCI.out.loci_merged_fa.collect()
    loci_merged_context_gz = MERGE_SEQS_LOCI.out.loci_merged_context_gz.collect()
    all_context_coords_bed = MERGE_SEQS_LOCI.out.all_context_coords_bed.collect()

	// Run reciprocal DIAMOND

		if (params.custom_reciprocal) {

			def reciprocal_ch = Channel.fromPath("${params.data_path}/${params.custom_reciprocal_db}", checkIfExists: true)

				// Run a file type check

					def reciprocal_type_ch = reciprocal_ch.map { file ->
						def extension = file.name.tokenize('.')[-1].toLowerCase()
						def type = extension in ['fasta', 'fa', 'fna'] ? 'fasta' : (extension == 'dmnd' ? 'dmnd' : null)
						if (type == null) {
							error ("Unsupported file extension: ${extension}")
						}
						return tuple(file, type)
					}

				// Branch the input channel based on file type

					reciprocal_type_ch.branch {
						fasta: it[1] == 'fasta'
						dmnd: it[1] == 'dmnd'
					}.set { reciprocal_branched }

				// Create reciprocal database if input fasta, else use provided dmnd db

					reciprocal_db = reciprocal_branched.fasta.map { file, type ->
						BUILD_RECIPROCAL("reciprocal", file)
					}.mix(reciprocal_branched.dmnd.map { file, type -> file })

        // run reciprocal DIAMOND and publish results
        SINGLE_RECIPROCAL_DIAMOND(reciprocal_db, loci_merged_fa)
        reciprocal_matches = SINGLE_RECIPROCAL_DIAMOND.out.reciprocal_matches.collect()
        reciprocal_seqs = SINGLE_RECIPROCAL_DIAMOND.out.reciprocal_seqs.collect()
        reciprocal_hits = SINGLE_RECIPROCAL_DIAMOND.out.reciprocal_hits.collect()

        // Find best hits    
        FIND_BEST_DIAMOND_HITS(forward_matches,
                        query_ch, reciprocal_matches,
                        reciprocal_seqs, reciprocal_hits)

        best_pairs_subsets = FIND_BEST_DIAMOND_HITS.out.best_pairs_txt.splitText(
                                                            by: params.pairs_per_task,
                                                            file: true)

        best_hit_proteins_val = FIND_BEST_DIAMOND_HITS.out.best_hits_fa_ch.collect()
        all_diamond_hits = FIND_BEST_DIAMOND_HITS.out.forward_plus_reciprocal_dmnd_hits.collect()
        
        // Make taxonomy and publish table for proteins
        hits_taxonomy = FETCH_HITS_TAXONOMY_FROM_ACCNS(all_diamond_hits)
    }

    else {
        // Reciprocal DIAMOND with rvdb and nr protein databases
        def reciprocal_nr_db_ch = Channel.fromPath("${params.data_path}/${params.reciprocal_nr_db}", checkIfExists: true)
        def reciprocal_rvdb_db_ch = Channel.fromPath("${params.data_path}/${params.reciprocal_rvdb_db}", checkIfExists: true)
        
        FULL_RECIPROCAL_DIAMOND(reciprocal_nr_db_ch, reciprocal_rvdb_db_ch,
                                loci_merged_fa, forward_matches,
                                query_ch)
        best_pairs_subsets = FULL_RECIPROCAL_DIAMOND.out.best_pairs_txt.splitText(
                                                            by: params.pairs_per_task,
                                                            file: true)
        best_hit_proteins_val = FULL_RECIPROCAL_DIAMOND.out.best_hits_fa_ch.collect()
        all_diamond_hits = FULL_RECIPROCAL_DIAMOND.out.mixed_hits.collect()

        // Read taxonomy table to build hits taxonomy
        def ncbi_tax_table_hits = Channel.fromPath("${params.data_path}/${params.ncbi_taxonomy_table}", checkIfExists: true)
        
        // Build hits taxonomy from annotated diamond database
        hits_taxonomy = BUILD_HITS_TAXONOMY_TABLE( FULL_RECIPROCAL_DIAMOND.out.reciprocal_nr_matches_ch,
                                                   FULL_RECIPROCAL_DIAMOND.out.reciprocal_rvdb_matches_ch,
                                                   ncbi_tax_table_hits)
        
    }
    
    // Reconstruction of encoded proteins. Returns number of STOP codons, frameshifts and indels
    GENEWISE(best_pairs_subsets,
            best_hit_proteins_val,
            loci_merged_fa,
            loci_merged_context_gz,
            all_context_coords_bed).collectFile(name: 'genewise.tsv',
                                            newLine: false,
                                            storeDir: "${params.outdir}/sql")

    // Publish files
    cat_forward = PUBLISH_FORWARD_DIAMOND(forward_matches, "forward-matches.dmnd.annot.tsv")
    predicted_orfs = ORF_EXTRACT(extract_seqs_outputs.context_fa_ch,
                                 extract_seqs_outputs.strict_coords_ch)

    cat_orfs = PUBLISH_PREDICTED_ORFS(predicted_orfs.orfs.collect(), "predicted_orfs.tsv")
    locus_assembly_maps = extract_seqs_outputs.locus_assembly_map_ch.collect()
    cat_assembly_map = PUBLISH_ASSEMBLY_MAP(locus_assembly_maps, "locus_assembly_map.tsv")

}
