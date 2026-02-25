/* 
----------------------------------------------------------------------------------------
Verify base dependencies & check required inputs
----------------------------------------------------------------------------------------
*/

// Verification workflow execution

workflow VERIFY {

	// Require nextflow

		if (!nextflow.version.matches('>=23.04.1')) {
			error ("ERROR: This workflow requires Nextflow version '23.04.1' or later. You are running '${nextflow.version}'.")
		}

	// Check assembly input mode and required parameters

		def supported_modes = ['ftp', 'local']
		def assembly_mode = (params.assembly_mode ?: 'ftp').toString().toLowerCase()

		if (!supported_modes.contains(assembly_mode)) {
			error "ERROR: Unsupported '--assembly_mode ${params.assembly_mode}'. Supported values: ftp, local."
		}

		if (assembly_mode == 'ftp' && !params.ftp_file) {
			error "ERROR: '--ftp_file' is required when '--assembly_mode ftp'."
		}

		if (assembly_mode == 'local' && !params.assembly_file) {
			error "ERROR: '--assembly_file' is required when '--assembly_mode local'."
		}

	// Check email only when Entrez-dependent features are enabled

		def needs_email = params.custom_reciprocal || assembly_mode == 'ftp' || params.get_all_metadata
		if (needs_email && !params.email) {
			error "ERROR: The '--email' parameter is required for Entrez API transactions in this configuration."
		}

	// Check FASTA inputs

		def fasta_extensions = ['fa', 'fna', 'fasta']
		def query_fasta = file("${params.data_path}/${params.query_file_aa}")

		if (query_fasta.exists()) {
			def file_extension = query_fasta.extension
				if (!fasta_extensions.contains(file_extension)) {
					error ("ERROR: File '${query_fasta}' does not have a '.fa', '.fna', or '.fasta' extension, please provide a FASTA file.")
				}
		} else {
			error ("ERROR: Reference file not found at '${query_fasta}'. Please check you've set the parameters '--data_path' and '--fasta' correctly.")
		}

	// Check local assembly file extension when local mode is used

		if (assembly_mode == 'local') {
			def local_assemblies = Channel.fromPath("${params.data_path}/${params.assembly_file}", checkIfExists: true).toList().getVal()
			if (local_assemblies.isEmpty()) {
				error ("ERROR: No assembly files matched '${params.data_path}/${params.assembly_file}'.")
			}

			def assembly_extensions = ['fa', 'fna', 'fasta', 'gz']
			local_assemblies.each { asm ->
				def ext = asm.extension?.toLowerCase()
				if (!assembly_extensions.contains(ext)) {
					error ("ERROR: Local assembly file '${asm}' must end with .fa, .fna, .fasta, or .gz.")
				}
			}
		}

}
