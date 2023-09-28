// Workflow specific params

// Import modules

// Run cloud workflow

workflow CLOUD {

    // Check that outdir does not exist. If not, proceed. If yes, warn.
    if (!file("$params.outdir").exists()) {

        // Run workflow
        println "CLOUD EXECUTION SCRIPT"

    }

    else {

        println ("Folder '$params.outdir' already exists. Remove or rename it before rerunning, or use --outdir to direct workflow outputs to an alternative directory.")

    }

}
