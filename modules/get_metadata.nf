process GET_METADATA {

	tag "${assembly_stats}"
	container 'oras://community.wave.seqera.io/library/biopython:1.85--a6f54362a38870a0'
	conda 'conda-forge::biopython=1.85'

	input:
	path assembly_stats

	output:
	path "assembly_metadata.tsv", emit: assembly_metadata
	publishDir "${params.outdir}/sql", mode: "copy", pattern: "assembly_metadata.tsv"

	script:
	"""

	get_assemblies_metadata.py ${assembly_stats} ${params.email} --outfile assembly_metadata.tsv
	
	"""

	stub:
	"""
	printf "hostName\tassembly_id\nunknown_host\tstub\n" > assembly_metadata.tsv
	"""

}
