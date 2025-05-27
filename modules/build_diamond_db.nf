process BUILD_DIAMOND_DB {

	tag "${sequences}"
	container 'oras://community.wave.seqera.io/library/diamond_seqkit_seqtk:6fc81cc10da8e7e4'
	conda 'bioconda::diamond=2.1.11'

	input:
	val label // label for output
	path sequences //input file

	output:
	path "${label}_db.dmnd", emit: vir_db_ch
	publishDir "${params.outdir}/${label}_db", mode: "copy"

	"""

	diamond makedb --in ${sequences} -d ${label}_db

	"""

}
