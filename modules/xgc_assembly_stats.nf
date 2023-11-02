process assembly_stats {

    input:
    path assembly

    output:
    stdout

    """

    stats.sh in=$assembly format=3 addname= | \
    grep -v n_scaffolds | \
    sed 's/\\/.*\\///g; s/_genomic.fna.gz//'

    """

}
