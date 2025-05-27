process FETCH_HOST_TAXONOMY {
	
	tag "${assembly_metadata}"
	container 'oras://community.wave.seqera.io/library/biopython_pandas:30b982b2baa37f61'
	conda 'conda-forge::biopython=1.85 conda-forge::pandas=2.2.3'

	input: 
	path assembly_metadata

	output:
	path "host_tax_information.tsv"
	publishDir "${params.outdir}/sql", mode: "copy", pattern: "host_tax_information.tsv"

	"""

	get_lineage_from_assembly_id.py ${assembly_metadata} ${params.email} --outfile host_tax_information.tsv
	
	"""

}
