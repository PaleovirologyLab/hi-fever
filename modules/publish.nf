process publish {

    input:
    path locus_assembly_map_collected
    path orfs_collected

    output:
    path "locus_assembly_map.txt"
    path "predicted_ORFs.txt"
    publishDir "${params.outdir}/sql", mode: "move", pattern: "locus_assembly_map.txt"
    publishDir "${params.outdir}/sql", mode: "move", pattern: "predicted_ORFs.txt"

    """

    cat $locus_assembly_map_collected > locus_assembly_map.txt
    cat $orfs_collected > predicted_ORFs.txt

    """

}
