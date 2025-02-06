/* 
----------------------------------------------------------------------------------------
Verify base dependencies & check required inputs
----------------------------------------------------------------------------------------
*/

// Verification workflow execution

workflow VERIFY {

	// Check required params are supplied before main workflow starts
		if (!params.email) {
			error "ERROR: The '--email' parameter is required for Entrez API transactions. Please provide a valid email address."
		}

}
