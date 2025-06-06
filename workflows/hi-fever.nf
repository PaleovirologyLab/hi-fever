// Parse inputs
include { PARSE_FTP } from '../modules/parse_ftp.nf'
include { CLUSTER_SEQS } from '../modules/cluster_seqs.nf'

// Process host genomes and get their taxonomical information
include { DOWNLOAD_ASSEMBLIES } from '../modules/download_assemblies.nf'
include { DOWNLOAD_EXTRACT_HOST_METADATA } from '../modules/process_host_info.nf'
include { GET_METADATA } from '../modules/get_metadata.nf'
include { ASSEMBLY_STATS } from '../modules/assembly_stats.nf'
include { FETCH_HOST_TAXONOMY } from '../modules/fetch_host_taxonomy.nf'
include { BUILD_HOST_TAXONOMY_TABLE } from '../modules/process_host_info.nf'
include { BLASTDB } from '../modules/blast_db.nf'

// DIAMOND-related process
include { BUILD_DIAMOND_DB as BUILD_QUERY} from '../modules/build_diamond_db.nf'
include { BUILD_DIAMOND_DB as BUILD_RECIPROCAL} from '../modules/build_diamond_db.nf'
include { FORWARD_DIAMOND } from '../modules/forward_diamond.nf'
include { SINGLE_RECIPROCAL_DIAMOND } from '../modules/diamond.nf'
include { FULL_RECIPROCAL_DIAMOND } from '../modules/diamond.nf'
include { FIND_BEST_DIAMOND_HITS } from '../modules/diamond.nf'

// Concatenate and extract best hits after forward DIAMOND
include { EXTRACT_SEQS_ANNOTATE_MATCHES } from '../modules/intersect_domains_merge_extract.nf'
include { MERGE_SEQS_LOCI} from '../modules/merge_seqs_loci.nf'

// Reconstruct interrupted EVE protein sequence
include { GENEWISE } from '../modules/genewise.nf'

// Extract orfs from forward DIAMOND
include { ORF_EXTRACT } from '../modules/orf_extract.nf'

// Getting hits taxonomy tables
include { BUILD_HITS_TAXONOMY_TABLE } from '../modules/hits_taxonomy.nf'
include { FETCH_HITS_TAXONOMY_FROM_ACCNS } from '../modules/hits_taxonomy.nf'

// Create summary table and extract fasta sequences for cdna and reconstructed proteins
include {CREATE_SUMMARY_TABLE_CUSTOM} from '../modules/create_summary_table.nf'
include {CREATE_SUMMARY_TABLE_FULL} from '../modules/create_summary_table.nf'

// Publish concatenated tables with results of all hosts
include { CONCATENATE_PUBLISH_TABLES as PUBLISH_PREDICTED_ORFS} from '../modules/utils.nf'
include { CONCATENATE_PUBLISH_TABLES as PUBLISH_FORWARD_DIAMOND} from '../modules/utils.nf'
include { CONCATENATE_PUBLISH_TABLES as PUBLISH_ASSEMBLY_MAP} from '../modules/utils.nf'

// Run workflow

workflow HIFEVER {

	// Create channels for protein FASTA query, assembly list (ftp links)

		def query_ch = Channel.fromPath("${params.data_path}/${params.query_file_aa}", checkIfExists: true)
		def ftp_ch = Channel.fromPath("${params.data_path}/${params.ftp_file}", checkIfExists: true)

	// If user provides own DMND query db (--query_db), create DIAMOND query channel from path

		if (params.query_db) {
			user_dmnd_db = Channel.fromPath("${params.data_path}/${params.user_dmnd_db}", checkIfExists: true)
		}

	// If clustering FASTA query, set query_ch to the output of clustering, else use input FASTA

		query_ch = params.cluster_query ? CLUSTER_SEQS(query_ch) : query_ch

	// DIAMOND query database: if user is providing it (--query_db), then use user_dmnd_db, else build new database from query_ch

		vir_db_ch = params.query_db ? user_dmnd_db : BUILD_QUERY("query", query_ch)

	// Unpack ftp list, download assemblies

		fetched_assembly_files = PARSE_FTP(ftp_ch) | flatten | DOWNLOAD_ASSEMBLIES

	// Add assembly accession as a meta field alongside assembly file path

		assembly_with_accession = fetched_assembly_files.map { assembly ->
										def fileName = assembly.baseName
										def accession = fileName.split('_')[0..1].join('_')
										def meta = [
										id: accession
										]
										return [meta, assembly]
								}

	// Collect statistics for downloaded assembly files

		assembly_stats = ASSEMBLY_STATS(fetched_assembly_files)
			.collectFile(name: 'assembly_stats.tsv', newLine: false, storeDir: "${params.outdir}/sql")

	// Get entrez metadata entries for either the downloaded assemblies or all eukaryotes

		if (!params.get_all_metadata) {
			// Download metadata only for genomes on ftp file
			metadata_channel = GET_METADATA(assembly_stats)
			FETCH_HOST_TAXONOMY(metadata_channel.assembly_metadata)

		} else {
			// Download assembly metadata for all eukaryotes
			DOWNLOAD_EXTRACT_HOST_METADATA()
			def ncbi_tax_table = Channel.fromPath("${params.data_path}/${params.ncbi_taxonomy_table}", checkIfExists: true)
			BUILD_HOST_TAXONOMY_TABLE( ftp_ch,
									DOWNLOAD_EXTRACT_HOST_METADATA.out.assembly_metadata_ch,
									ncbi_tax_table)
		}

	// Run forward DIAMOND using chunks of the genome as queries against the viral DMND database

		forward_diamond_out = FORWARD_DIAMOND(assembly_with_accession.combine(vir_db_ch))

	// Make BLAST db for assembly

		BLASTDB(assembly_with_accession)

	// Rejoin dmnd.nsq with nsq db

		rejoinDb = forward_diamond_out.join(BLASTDB.out, by: [0]) // Joins on shared metadata (assembly)

	// Extract genome locus with hits, and their flanking regions

		extract_seqs_outputs = EXTRACT_SEQS_ANNOTATE_MATCHES(rejoinDb)

	// Inputs for reciprocal DIAMOND

		forward_matches = extract_seqs_outputs.forward_matches.collect()
		strict_fastas_collected = extract_seqs_outputs.strict_fa_ch.collect()
		context_fastas_collected = extract_seqs_outputs.context_fa_ch.collect()

	// Concatenate sequences into single file

		MERGE_SEQS_LOCI(strict_fastas_collected, context_fastas_collected)
		loci_merged_fa = MERGE_SEQS_LOCI.out.loci_merged_fa.collect()
		loci_merged_context_gz = MERGE_SEQS_LOCI.out.loci_merged_context_gz.collect()
		all_context_coords_bed = MERGE_SEQS_LOCI.out.all_context_coords_bed.collect()


	// Run reciprocal DIAMOND

		// Workflow for user custom reciprocal database	(single db)

			if (params.custom_reciprocal) {

				// Create input channel, check type

				reciprocal_ch = Channel.fromPath("${params.data_path}/${params.custom_reciprocal_db}", checkIfExists: true)
					.map { file ->
						def extension = file.extension.toLowerCase()
						def type
						if (extension == 'fa' || extension == 'fasta' || extension == 'fna') {
							type = 'fasta'
						} else if (extension == 'dmnd') {
							type = 'dmnd'
						} else {
							error "Unsupported database file extension: ${extension}. Supported formats: .fa, .fasta, .fna, .dmnd"
						}
						return [ file, type ]
					}

				// Branch channel based on file type & build reciprocal database if FASTA

				reciprocal_ch.branch {
					fasta: it[1] == 'fasta'
					dmnd: it[1] == 'dmnd'
				}.set { reciprocal_branched_ch }

				// Set reciprocal db channel - if dmdn, BUILD_RECIPROCAL will not start, if fasta, it will build reciprocal database

				reciprocal_db = reciprocal_branched_ch.dmnd.map { file, type -> file }.mix(BUILD_RECIPROCAL("reciprocal", reciprocal_branched_ch.fasta.map { file, type -> file}))

				// run reciprocal DIAMOND and publish results

				SINGLE_RECIPROCAL_DIAMOND(reciprocal_db, loci_merged_fa)
				reciprocal_matches = SINGLE_RECIPROCAL_DIAMOND.out.reciprocal_matches
				reciprocal_seqs = SINGLE_RECIPROCAL_DIAMOND.out.reciprocal_seqs
				reciprocal_hits = SINGLE_RECIPROCAL_DIAMOND.out.best_reciprocal_hits
                all_reciprocal_hits = SINGLE_RECIPROCAL_DIAMOND.out.reciprocal_hits

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
				hits_taxonomy = FETCH_HITS_TAXONOMY_FROM_ACCNS(all_reciprocal_hits)

		} else {

		// Workflow for reciprocal DIAMOND with rvdb and nr protein databases (i.e. full)

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
			hits_taxonomy = BUILD_HITS_TAXONOMY_TABLE(FULL_RECIPROCAL_DIAMOND.out.reciprocal_nr_matches_ch,
													FULL_RECIPROCAL_DIAMOND.out.reciprocal_rvdb_matches_ch,
													ncbi_tax_table_hits)

		}

	// Reconstruction of encoded proteins. Returns number of STOP codons, frameshifts and indels

        genewise_result = GENEWISE(best_pairs_subsets,
                best_hit_proteins_val,
                loci_merged_fa,
                loci_merged_context_gz,
                all_context_coords_bed)


    // Merge genewise results into a single file            
       
        merged_genewise_file = genewise_result.genewise_file.collectFile(name: 'genewise.tsv', 
                                                                        newLine: false, 
                                                                        storeDir: "${params.outdir}/sql")

	// Publish files

		cat_forward = PUBLISH_FORWARD_DIAMOND(forward_matches, "forward-matches.dmnd.annot.tsv")
		predicted_orfs = ORF_EXTRACT(extract_seqs_outputs.context_fa_ch,
									extract_seqs_outputs.strict_coords_ch)

		cat_orfs = PUBLISH_PREDICTED_ORFS(predicted_orfs.orfs.collect(), "predicted_orfs.tsv")
		locus_assembly_maps = extract_seqs_outputs.locus_assembly_map_ch.collect()
		cat_assembly_map = PUBLISH_ASSEMBLY_MAP(locus_assembly_maps, "locus_assembly_map.tsv")

    
    // Get a summary table with statistics and classification of query loci
        if (params.custom_reciprocal) {
                summary_table = CREATE_SUMMARY_TABLE_CUSTOM(
                    reciprocal_rvdb = reciprocal_matches.collect(),
                    taxonomy = hits_taxonomy,
                    assembly_map = cat_assembly_map,
                    assembly_metadata = metadata_channel.assembly_metadata,
                    genewise = merged_genewise_file
            )
        } else {
            summary_table = CREATE_SUMMARY_TABLE_FULL(
                reciprocal_nr = FULL_RECIPROCAL_DIAMOND.out.reciprocal_nr_matches_ch,
                reciprocal_rvdb = FULL_RECIPROCAL_DIAMOND.out.reciprocal_rvdb_matches_ch,
                taxonomy = hits_taxonomy,
                assembly_map = cat_assembly_map,
                assembly_metadata = metadata_channel.assembly_metadata,
                genewise = merged_genewise_file
            )

        }


}
