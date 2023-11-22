process get_assembly_metadata {

    output:
    path "assembly_metadata.tsv"
    publishDir "${params.outdir}/sql", mode: "move", pattern: "assembly_metadata.tsv"

    """

    # Download current eukaryotic assembly summaries

    wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/assembly_summary.txt; mv assembly_summary.txt rs_assembly_summary_invert.txt
    wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt; mv assembly_summary.txt rs_assembly_summary_fungi.txt
    wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/assembly_summary.txt; mv assembly_summary.txt rs_assembly_summary_plant.txt
    wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/assembly_summary.txt; mv assembly_summary.txt rs_assembly_summary_protozoa.txt
    wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/assembly_summary.txt; mv assembly_summary.txt rs_assembly_summary_mammal.txt
    wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/assembly_summary.txt; mv assembly_summary.txt rs_assembly_summary_vert_other.txt
    wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/invertebrate/assembly_summary.txt; mv assembly_summary.txt gb_assembly_summary_invert.txt
    wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/assembly_summary.txt; mv assembly_summary.txt gb_assembly_summary_fungi.txt
    wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/assembly_summary.txt; mv assembly_summary.txt gb_assembly_summary_plant.txt
    wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/protozoa/assembly_summary.txt; mv assembly_summary.txt gb_assembly_summary_protozoa.txt
    wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/assembly_summary.txt; mv assembly_summary.txt gb_assembly_summary_mammal.txt
    wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_other/assembly_summary.txt; mv assembly_summary.txt gb_assembly_summary_vert_other.txt

    # Remove header lines and concatenate

    for i in ??_assembly_summary_*; do grep -v "README_assembly_summary\\|#assembly_accession" \$i >> metadata_partial.txt; done

    # Extract AssemblyID

    cut -f20 metadata_partial.txt | sed 's/.*\\///g' > assemblyID

    # Paste

    paste assemblyID metadata_partial.txt > assembly_metadata.tsv

    """

}