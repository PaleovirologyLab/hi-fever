process BLASTDB {

	// Directives

	debug false
	container 'oras://community.wave.seqera.io/library/blast:2.16.0--ee6ae29ad5529d04'
	conda 'bioconda::blast=2.16.0'
	tag "${meta.id}"

	input:
	tuple val(meta), path(assembly)

	output:
	tuple val(meta), path("*.nsq")

	"""

    if [[ "${assembly}" == *.gz ]]; then
        gunzip -c ${assembly}
    else
        cat ${assembly}
    fi | makeblastdb -in - -out ${meta.id} -title ${meta.id} -dbtype nucl -parse_seqids

	"""

}
