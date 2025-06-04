process BLASTDB {

	// Directives

	debug false
	container 'oras://community.wave.seqera.io/library/blast:2.16.0--ee6ae29ad5529d04'
	conda 'bioconda::blast=2.16.0'
	tag "${meta.id}"

	input:
	tuple val(meta), path(assembly)

	output:
	tuple val(meta), path("*.gz*nsq")

	"""

    gunzip -c ${assembly} | makeblastdb -in - -out ${assembly} -title ${assembly} -dbtype nucl -parse_seqids

	"""

}
