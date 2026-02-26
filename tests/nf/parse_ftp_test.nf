nextflow.enable.dsl=2

include { PARSE_FTP } from '../../modules/parse_ftp.nf'

workflow {
    ftp_ch = Channel.fromPath(params.ftp_input, checkIfExists: true)
    parsed = PARSE_FTP(ftp_ch)

    parsed
        .flatten()
        .map { it.getName() }
        .collectFile(name: 'parsed_manifest.txt', newLine: true, storeDir: params.outdir)
}
