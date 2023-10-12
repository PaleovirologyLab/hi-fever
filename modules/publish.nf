process publish {

    input:
    path strict_fastas_collected
    path context_fastas_collected
    path locus_assembly_map_collected
    path orfs_collected
    path annotated_hits_collected

    output:
    path "loci-merged-coordinates.fasta.gz"
    path "loci-context-coordinates.fasta.gz"
    path "locus_assembly_map.txt"
    path "predicted_ORFs.txt"
    path "matches.dmnd.annot.tsv"
    publishDir "${params.outdir}/accessory_fastas", mode: "move", pattern: "*.fasta.gz"
    publishDir "${params.outdir}/sql", mode: "move", pattern: "locus_assembly_map.txt"
    publishDir "${params.outdir}/sql", mode: "move", pattern: "predicted_ORFs.txt"
    publishDir "${params.outdir}/sql", mode: "move", pattern: "matches.dmnd.annot.tsv"

    """

    cat $strict_fastas_collected | gzip > loci-merged-coordinates.fasta.gz
    cat $context_fastas_collected | gzip > loci-context-coordinates.fasta.gz
    cat $locus_assembly_map_collected > locus_assembly_map.txt
    cat $orfs_collected > predicted_ORFs.txt
    cat $annotated_hits_collected > matches.dmnd.annot.tsv

    """

}
