#!/home/user/nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hi-fever (high-throughput nextflow EVE recovery)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GitHub: https://github.com/Paleovirology/hi-fever

Maintainers:
- Cormac Kinsella (cormac.kinsella@evobio.eu)
- José Gabriel Niño Barreat (jose.ninobarreat@biology.ox.ac.uk)
- Emma Harding (emma.harding@biology.ox.ac.uk)
- Laura Munoz-Buena 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Syntax version

nextflow.enable.dsl=2

// Import workflow scripts

include { HIFEVER } from './workflows/hi-fever.nf'

//  Workflow definition
workflow {

    HIFEVER ()
	
}