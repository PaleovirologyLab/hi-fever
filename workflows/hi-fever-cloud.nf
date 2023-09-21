// Import modules


// Run cloud workflow

workflow CLOUD {

    // Check that outdir does not exist. If not, proceed. If yes, warn.
    if (!file("$params.outdir").exists()) {

        // Run workflow
        println "CLOUD EXECUTION SCRIPT"

    }

    else {

        println ("Folder '$params.outdir' exists. Either remove it before rerunning, or rename to preserve the contents, e.g., YYYY-MM-DD-experiment_description")

    }

}
