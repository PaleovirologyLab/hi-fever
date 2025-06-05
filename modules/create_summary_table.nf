process CREATE_SUMMARY_TABLE{
        // container 'oras://community.wave.seqera.io/library/biopython_modin:2303fabbcd10e9be'
	    // conda 'conda-forge::modin=0.32.0 conda-forge::biopython=1.85'


    input:
        path reciprocal_nr
        path reciprocal_rvdb
        path taxonomy
        path assembly_map
        path assembly_metadata
        path genewise

    output:
        path "*_classified.tsv", emit: summary_table
        path "*_aa_genewise.fasta", emit: aa_fasta
        path "*_cdna_genewise.fasta", emit: cdna_fasta

        publishDir "${params.outdir}/sql", mode: "copy", pattern: "*.tsv"
        publishDir "${params.outdir}/sql", mode: "copy", pattern: "*.fasta"

        """
        create_summary_table.py \\
        --reciprocal_nr ${reciprocal_nr} \\
        --reciprocal_rvdb ${reciprocal_rvdb} \\
        --taxonomy ${taxonomy} \\
        --assembly_map ${assembly_map} \\
        --assembly_metadata ${assembly_metadata} \\
        --genewise ${genewise} \\
        
        """

}