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

	// Check required params are supplied before main workflow starts

		if (!params.email) {
			error "ERROR: The '--email' parameter is required for Entrez API transactions. Please provide a valid email address."
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

}
