

process SINGLE_RECIPROCAL_DIAMOND {

	container 'oras://community.wave.seqera.io/library/diamond_seqkit_seqtk:6fc81cc10da8e7e4'
	conda 'bioconda::diamond=2.1.11'


    input:
    path reciprocal_db
    path loci_merged_fa

    output:
    path "reciprocal-matches.dmnd.tsv", emit: reciprocal_matches
    path "reciprocal_hits.txt", emit: reciprocal_hits
    path "reciprocal_seqs.fasta", emit: reciprocal_seqs

    publishDir "${params.outdir}/sql", mode: "copy", pattern: "reciprocal-matches.dmnd.tsv"
    

    """ 

	# Run reciprocal nr search, extract best reciprocal hit pairs & their alignment length

    diamond blastx \
    --${params.diamond_mode} \
    --matrix ${params.diamond_matrix} \
    --masking seg \
    -d ${reciprocal_db} \
    -q ${loci_merged_fa} \
    -e 1e-5 \
    -k 10 \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle full_sseq | \
    tee >(sort -k1,1 -k12,12nr | sort -u -k1,1 | \
    tee >(awk 'BEGIN{OFS="\t"}; {print \$2, \$1, \$4 * 3, "reciprocal"}' > reciprocal_hits.txt) \
    | cut -f2,19 | sort -u -k1,1 | awk '{print ">"\$1"\\n"\$2}' > reciprocal_seqs.fasta) | \
    cut -f1-18 > \
    reciprocal-matches.dmnd.tsv

    # Zip up query file for outdir publication
    cat ${loci_merged_fa} > loci-merged-coordinates.fasta.temp
    gzip loci-merged-coordinates.fasta.temp > loci-merged-coordinates.fasta.gz
    """

}

process FIND_BEST_DIAMOND_HITS {

	container 'oras://community.wave.seqera.io/library/diamond_seqkit_seqtk:6fc81cc10da8e7e4'
	conda 'bioconda::seqtk=r93'

    input:
    path forward_matches
    path query_proteins
    path reciprocal_matches
    path reciprocal_seqs
    path reciprocal_hits

    output:
    path "best_hits.fasta", emit: best_hits_fa_ch
    path "best_pairs.txt", emit: best_pairs_txt
    path "mixed_hits.txt", emit: forward_plus_reciprocal_dmnd_hits
    
    // publishDir "${params.outdir}", mode: "copy", pattern: "best_hits.fasta"
    // publishDir "${params.outdir}", mode: "copy", pattern: "best_pairs.txt"
    // publishDir "${params.outdir}", mode: "copy", pattern: "mixed_hits.txt"
    // publishDir "${params.outdir}/sql", mode: "copy", pattern: "all_forward.dmnd.annot.tsv"

    """
    cat ${forward_matches} > all_forward.dmnd.annot.tsv

    # Extract best forward hit pairs & their alignment length
    awk 'BEGIN{OFS="\t"}; {print \$1,\$4,\$9-\$8, "forward"}' all_forward.dmnd.annot.tsv | \
    uniq > \
    forward_hits.txt
    
    # Merge reciprocal and forward hits
    cat forward_hits.txt ${reciprocal_hits} > mixed_hits.txt

    # Determine best single protein-locus pair from both the forward & reciprocal searches

    sort -k2,2 -k3,3nr mixed_hits.txt | \
    sort -u -k2,2 | \
    tee >(grep reciprocal | cut -f1 | sort | uniq > extract_from_reciprocal.txt) | \
    tee >(grep forward | cut -f1 | sort | uniq > extract_from_forward.txt) > \
    best_pairs_partial.txt

    # Add true genomic start for later conversion

    cut -f2 best_pairs_partial.txt | sed 's/[:-]/\t/g' > coords.txt
    paste best_pairs_partial.txt coords.txt > best_pairs.txt

    # Extract best proteins

    seqtk subseq ${reciprocal_seqs} extract_from_reciprocal.txt | \
    sed 's/ .*//' > \
    best_hits.fasta

    seqtk subseq ${query_proteins} extract_from_forward.txt | \
    sed 's/ .*//' >> \
    best_hits.fasta

    """

}


process FULL_RECIPROCAL_DIAMOND {

	container 'oras://community.wave.seqera.io/library/diamond_seqkit_seqtk:6fc81cc10da8e7e4'
	conda 'bioconda::diamond=2.1.11 bioconda::seqtk=r93'

    input:
    path reciprocal_nr_db
    path reciprocal_rvdb_db
    path loci_merged_fa
    path forward_matches
    path query_proteins


    output:
    path "reciprocal-nr-matches.dmnd.tsv", emit: reciprocal_nr_matches_ch
    path "reciprocal-rvdb-matches.dmnd.tsv", emit: reciprocal_rvdb_matches_ch
    path "best_hits.fasta", emit: best_hits_fa_ch
    path "best_pairs.txt", emit: best_pairs_txt
    path "mixed_hits.txt", emit: mixed_hits
 
    publishDir "${params.outdir}/sql", mode: "copy", pattern: "reciprocal-nr-matches.dmnd.tsv"
    publishDir "${params.outdir}/sql", mode: "copy", pattern: "reciprocal-rvdb-matches.dmnd.tsv"
    
    """ 
    # Run reciprocal nr search, extract best reciprocal hit pairs & their alignment length
    diamond blastx \
    --${params.diamond_mode} \
    --matrix ${params.diamond_matrix} \
    --masking seg \
    -d ${reciprocal_nr_db} \
    -q ${loci_merged_fa} \
    -e 1e-5 \
    -k 10 \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms skingdoms sphylums stitle full_sseq | \
    tee >(sort -k1,1 -k12,12nr | sort -u -k1,1 | tee >(awk 'BEGIN{OFS="\t"}; {print \$2, \$1, \$4 * 3, "reciprocal-nr"}' > mixed_hits.txt) | cut -f2,19 | sort -u -k1,1 | awk '{print ">"\$1"\\n"\$2}' > reciprocal_nr_subset.fasta) | \
    cut -f1-18 > \
    reciprocal-nr-matches.dmnd.tsv

    # Run reciprocal RVDB search, extract best reciprocal hit pairs & their alignment length

    diamond blastx \
    --${params.diamond_mode} \
    --matrix ${params.diamond_matrix} \
    --masking seg \
    -d ${reciprocal_rvdb_db} \
    -q ${loci_merged_fa} \
    -e 1e-5 \
    -k 10 \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms skingdoms sphylums stitle full_sseq | \
    tee >(sort -k1,1 -k12,12nr | sort -u -k1,1 | tee >(awk 'BEGIN{OFS="\t"}; {print \$2, \$1, \$4 * 3, "reciprocal-rvdb"}' >> mixed_hits.txt) | cut -f2,19 | sort -u -k1,1 | awk '{print ">"\$1"\\n"\$2}' > reciprocal_rvdb_subset.fasta) | \
    cut -f1-18 > \
    reciprocal-rvdb-matches.dmnd.tsv

    # Concatenate forward hit tables

    cat ${forward_matches} > matches.dmnd.annot.tsv

    # Extract best forward hit pairs & their alignment length

    awk 'BEGIN{OFS="\t"}; {print \$1,\$4,\$9-\$8, "forward"}' matches.dmnd.annot.tsv | \
    uniq >> \
    mixed_hits.txt

    # Determine best single protein-locus pair from both the forward & reciprocal searches

    sort -k2,2 -k3,3nr mixed_hits.txt | \
    sort -u -k2,2 | \
    tee >(grep reciprocal-nr | cut -f1 | sort | uniq > extract_from_reciprocal-nr.txt) | \
    tee >(grep reciprocal-rvdb | cut -f1 | sort | uniq > extract_from_reciprocal-rvdb.txt) | \
    tee >(grep forward | cut -f1 | sort | uniq > extract_from_forward.txt) > \
    best_pairs_partial.txt

    # Add true genomic start for later conversion

    cut -f2 best_pairs_partial.txt | sed 's/[:-]/\t/g' > coords.txt
    paste best_pairs_partial.txt coords.txt > best_pairs.txt

    # Extract best proteins

    seqtk subseq reciprocal_nr_subset.fasta extract_from_reciprocal-nr.txt | \
    sed 's/ .*//' > \
    best_hits.fasta

    seqtk subseq reciprocal_rvdb_subset.fasta extract_from_reciprocal-rvdb.txt | \
    sed 's/ .*//' >> \
    best_hits.fasta

    seqtk subseq ${query_proteins} extract_from_forward.txt | \
    sed 's/ .*//' >> \
    best_hits.fasta

    """

}
