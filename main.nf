#!/home/user/nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hi-fever (high-throughput nextflow EVE recovery)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GitHub: https://github.com/Paleovirology/hi-fever

Maintainers:
- Cormac Kinsella (cormac.kinsella@evobio.eu)
- José Gabriel Niño Barreat (jose.ninobarreat@biology.ox.ac.uk)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Syntax version

nextflow.enable.dsl=2

// Import workflow scripts

include { HIFEVER } from './workflows/hi-fever.nf'
include { CLOUD } from './workflows/hi-fever-cloud.nf'

//  Workflow definition

workflow {
    if ( "$params.entry" == "LOCAL" ) {
    HIFEVER ()
    }

    else if ( "$params.entry" == "CLOUD") {
    CLOUD ()
    }
}