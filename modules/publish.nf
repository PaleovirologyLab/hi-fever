process publish {

    input:
    path strict_fastas_collected
    path context_fastas_collected
    path locus_assembly_map_collected

    output:
    path "loci-merged-coordinates.fasta.gz"
    path "loci-context-coordinates.fasta.gz"
    path "locus_assembly_map.txt"
    publishDir "${params.outdir}/accessory_fastas", mode: "move", pattern: "*.fasta.gz"
    publishDir "${params.outdir}", mode: "move", pattern: "locus_assembly_map.txt"

    """

    cat $strict_fastas_collected | gzip > loci-merged-coordinates.fasta.gz
    cat $context_fastas_collected | gzip > loci-context-coordinates.fasta.gz
    cat $locus_assembly_map_collected > locus_assembly_map.txt

    """

}
