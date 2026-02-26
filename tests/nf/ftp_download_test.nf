nextflow.enable.dsl=2

include { PARSE_FTP } from '../../modules/parse_ftp.nf'
include { DOWNLOAD_ASSEMBLIES } from '../../modules/download_assemblies.nf'

workflow {
    ftp_ch = Channel.fromPath(params.ftp_input, checkIfExists: true)
    parsed = PARSE_FTP(ftp_ch)
    downloads = DOWNLOAD_ASSEMBLIES(parsed.flatten())

    downloads
        .map { it.getName() }
        .collectFile(name: 'downloaded_manifest.txt', newLine: true, storeDir: params.outdir)
}
