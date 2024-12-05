process diamond {

    maxForks params.diamond_forks

    input:
    tuple path(assembly), path(db)

    output:
    tuple path("*.dmnd.tsv"), path("*.gz*nsq")
    publishDir "${params.outdir}", mode: "copy", pattern: "*.dmnd.tsv"
    publishDir "${params.outdir}", mode: "copy", pattern: "*.gz*nsq"

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
    > matches.reformatted.dmnd.tsv

    # rm \$chunks

    # rm "\$(readlink -f \$assembly)"

    """

}

process build_diamond_db {
    input:
    val label // label for output
    path sequences //input file

    output:
    //path "DB_clu_rep.fasta", emit: clust_ch
    path "${label}_db.dmnd", emit: vir_db_ch
    publishDir "${params.outdir}/${label}_db", mode: "copy"

    """

    diamond makedb --in $sequences -d ${label}_db

    """

}

process reciprocal_diamond_single {

    input:
    path reciprocal_db
    path loci_merged_fa
    // path strict_fastas_collected
    // path context_fastas_collected
    // path reciprocal_rvdb_db
    // path annotated_hits_collected
    // path query_proteins

    output:
    path "reciprocal-matches.dmnd.tsv", emit: reciprocal_matches
    path "reciprocal_hits.txt", emit: reciprocal_hits
    path "reciprocal_seqs.fasta", emit: reciprocal_seqs
    path "loci-merged-coordinates.fasta.gz", emit: loci_merged_fa_gz
    publishDir "${params.outdir}", mode: "copy", pattern: "reciprocal-matches.dmnd.tsv"
    publishDir "${params.outdir}", mode: "copy", pattern: "reciprocal_hits.txt"
    publishDir "${params.outdir}", mode: "copy", pattern: "reciprocal_seqs.fasta"
    publishDir "${params.outdir}", mode: "copy", pattern: "loci-merged-coordinates.fasta.gz"
    

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

process reciprocal_diamond_double {

    input:
    path reciprocal_nr_db
    path reciprocal_rvdb_db
    path loci_merged_fa
    // path strict_fastas_collected
    // path context_fastas_collected
    // path reciprocal_rvdb_db
    // path annotated_hits_collected
    // path query_proteins

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

    """

}