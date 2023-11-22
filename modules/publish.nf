process publish {

    input:
    path locus_assembly_map_collected
    path orfs_collected

    output:
    path "locus_assembly_map.tsv"
    path "predicted_ORFs.tsv"
    publishDir "${params.outdir}/sql", mode: "move", pattern: "locus_assembly_map.tsv"
    publishDir "${params.outdir}/sql", mode: "move", pattern: "predicted_ORFs.tsv"

    """

    cat $locus_assembly_map_collected > locus_assembly_map.tsv
    cat $orfs_collected > predicted_ORFs.tsv

    """

}
