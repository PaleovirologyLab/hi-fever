process ASSEMBLY_STATS {

	tag "${assembly}"
	container 'oras://community.wave.seqera.io/library/bbmap:39.15--65f46923e9870921'
	conda 'bioconda::bbmap=39.15'

	input:
	path assembly

	output:
	stdout

	"""

	stats.sh in=${assembly} format=3 addname= | \
	grep -v n_scaffolds | \
	sed 's/\\/.*\\///g; s/_genomic.fna.gz//'

	"""

}
