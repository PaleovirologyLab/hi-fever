process build_diamond_db {
    input:
    val label // label for output
    path sequences //input file

    output:
    path "${label}_db.dmnd", emit: vir_db_ch
    publishDir "${params.outdir}/${label}_db", mode: "copy"

    """
    diamond makedb --in $sequences -d ${label}_db

    """

}

process forward_diamond {

    maxForks params.diamond_forks

    input:
    tuple path(assembly), path(db)

    output:
    tuple path("*.dmnd.tsv"), path("*.gz*nsq")
    publishDir "${params.outdir}", mode: "copy", pattern: "*.dmnd.tsv"
    // publishDir "${params.outdir}", mode: "copy", pattern: "*.gz*nsq"

    """

    db=$db
    assembly=$assembly
    chunks=\$(echo \$assembly | sed 's/_genomic.*/_genomic_chunks.fna.gz/')
    cpu_count=\$(awk -v total_cpu=\$(nproc) 'BEGIN {printf "%.0f\\n", (total_cpu > 1) ? total_cpu / $params.diamond_forks : 1}')

    seqkit sliding -s $params.chunk_size -W $params.chunk_size -g $assembly -o \$chunks

    diamond blastx \
    --$params.diamond_mode \
    --matrix $params.diamond_matrix \
    --masking seg \
    -d \$db \
    -q \$chunks \
    -o matches.out \
    -p \$cpu_count \
    --max-target-seqs $params.diamond_max_target_seqs \
    --outfmt 6 qseqid qstart qend qframe qlen sseqid sstart send slen evalue bitscore pident length mismatch gapopen &

    gunzip -c \$assembly | makeblastdb -in - -out \${assembly/*\\//} -title \${assembly/*\\//} -dbtype nucl -parse_seqids &

    wait

    sed 's/_sliding:/\\t/' matches.out | sed 's/-/\\t/' | \
    awk -v OFS='\\t' '{print \$1,\$2+\$4-1,\$2+\$5-1,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15,\$16,\$17}' \
    > "\${assembly}_forward-matches-raw.dmnd.tsv"

    # rm \$chunks

    # rm "\$(readlink -f \$assembly)"

    """

}

process reciprocal_diamond {

    input:
    path reciprocal_db
    path loci_merged_fa

    output:
    path "reciprocal-matches.dmnd.tsv", emit: reciprocal_matches
    path "reciprocal_hits.txt", emit: reciprocal_hits
    path "reciprocal_seqs.fasta", emit: reciprocal_seqs

    publishDir "${params.outdir}", mode: "copy", pattern: "reciprocal-matches.dmnd.tsv"
    publishDir "${params.outdir}", mode: "copy", pattern: "reciprocal_hits.txt"
    publishDir "${params.outdir}", mode: "copy", pattern: "reciprocal_seqs.fasta"
    

    """ 
   # Run reciprocal nr search, extract best reciprocal hit pairs & their alignment length

    diamond blastx \
    --$params.diamond_mode \
    --matrix $params.diamond_matrix \
    --masking seg \
    -d $reciprocal_db \
    -q $loci_merged_fa \
    -e 1e-5 \
    -k 10 \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle full_sseq | \
    tee >(sort -k1,1 -k12,12nr | sort -u -k1,1 | \
    tee >(awk 'BEGIN{OFS="\t"}; {print \$2, \$1, \$4 * 3, "reciprocal"}' > reciprocal_hits.txt) \
    | cut -f2,19 | sort -u -k1,1 | awk '{print ">"\$1"\\n"\$2}' > reciprocal_seqs.fasta) | \
    cut -f1-18 > \
    reciprocal-matches.dmnd.tsv

    # Zip up query file for outdir publication
    cat $loci_merged_fa > loci-merged-coordinates.fasta.temp
    gzip loci-merged-coordinates.fasta.temp > loci-merged-coordinates.fasta.gz
    """

}

process find_best_diamond_hits {

    input:
    path forward_matches
    path query_proteins
    path reciprocal_matches
    path reciprocal_seqs
    path reciprocal_hits

    output:
    path "best_hits.fasta", emit: best_hits_fa_ch
    path "best_pairs.txt", emit: best_pairs_txt
    
    // publishDir "${params.outdir}", mode: "copy", pattern: "best_hits.fasta"
    // publishDir "${params.outdir}", mode: "copy", pattern: "best_pairs.txt"

    """
    cat $forward_matches > forward_matches.dmnd.annot.tsv

    # Extract best forward hit pairs & their alignment length
    awk 'BEGIN{OFS="\t"}; {print \$1,\$4,\$9-\$8, "forward"}' forward_matches.dmnd.annot.tsv | \
    uniq > \
    forward_hits.txt
    
    # Merge reciprocal and forward hits
    cat forward_hits.txt $reciprocal_hits > mixed_hits.txt

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

    seqtk subseq $reciprocal_seqs extract_from_reciprocal.txt | \
    sed 's/ .*//' > \
    best_hits.fasta

    seqtk subseq $query_proteins extract_from_forward.txt | \
    sed 's/ .*//' >> \
    best_hits.fasta

    """

}


process reciprocal_diamond_full {

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
 
    publishDir "${params.outdir}/sql", mode: "copy", pattern: "reciprocal-nr-matches.dmnd.tsv"
    publishDir "${params.outdir}/sql", mode: "copy", pattern: "reciprocal-rvdb-matches.dmnd.tsv"
    
    """ 
    # Run reciprocal nr search, extract best reciprocal hit pairs & their alignment length
    diamond blastx \
    --$params.diamond_mode \
    --matrix $params.diamond_matrix \
    --masking seg \
    -d $reciprocal_nr_db \
    -q $loci_merged_fa \
    -e 1e-5 \
    -k 10 \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms skingdoms sphylums stitle full_sseq | \
    tee >(sort -k1,1 -k12,12nr | sort -u -k1,1 | tee >(awk 'BEGIN{OFS="\t"}; {print \$2, \$1, \$4 * 3, "reciprocal-nr"}' > mixed_hits.txt) | cut -f2,19 | sort -u -k1,1 | awk '{print ">"\$1"\\n"\$2}' > reciprocal_nr_subset.fasta) | \
    cut -f1-18 > \
    reciprocal-nr-matches.dmnd.tsv

    # Run reciprocal RVDB search, extract best reciprocal hit pairs & their alignment length

    diamond blastx \
    --$params.diamond_mode \
    --matrix $params.diamond_matrix \
    --masking seg \
    -d $reciprocal_rvdb_db \
    -q $loci_merged_fa \
    -e 1e-5 \
    -k 10 \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms skingdoms sphylums stitle full_sseq | \
    tee >(sort -k1,1 -k12,12nr | sort -u -k1,1 | tee >(awk 'BEGIN{OFS="\t"}; {print \$2, \$1, \$4 * 3, "reciprocal-rvdb"}' >> mixed_hits.txt) | cut -f2,19 | sort -u -k1,1 | awk '{print ">"\$1"\\n"\$2}' > reciprocal_rvdb_subset.fasta) | \
    cut -f1-18 > \
    reciprocal-rvdb-matches.dmnd.tsv

    # Concatenate forward hit tables

    cat $forward_matches > matches.dmnd.annot.tsv

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

    seqtk subseq $query_proteins extract_from_forward.txt | \
    sed 's/ .*//' >> \
    best_hits.fasta

    """

}
