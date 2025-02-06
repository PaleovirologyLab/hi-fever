#!/home/user/nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hi-fever (high-throughput nextflow EVE recovery)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GitHub: https://github.com/PaleovirologyLab/hi-fever

Maintainers:
- Cormac Kinsella (cormac.kinsella@scilifelab.se)
- José Gabriel Niño Barreat (jose.ninobarreat@biology.ox.ac.uk)
- Emma Harding (emma.harding@biology.ox.ac.uk)
- Laura Munoz-Buena 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Import workflow scripts

include { VERIFY } from './subworkflows/verify.nf'
include { HIFEVER } from './workflows/hi-fever.nf'

//  Workflow definition
workflow {

    VERIFY ()
    HIFEVER ()
	
}