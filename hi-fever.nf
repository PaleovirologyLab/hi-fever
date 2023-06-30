#!/home/user/nextflow

// Syntax version

nextflow.enable.dsl=2

// Script parameters

params.ftp_file = "$PWD/ftp_list.txt"
params.query_file_aa = "$PWD/protein_query.fasta"
params.phmms = "$PWD/domains-v*"
params.mmseqs_minseqid = "0.95"
params.mmseqs_cover = "0.90"
params.diamond_mode = "very-sensitive"
params.diamond_matrix = "BLOSUM62"
params.diamond_cpus = "12"

// Define workflow processes

process build_db {

    output:
    path "DB_clu_rep.fasta", emit: clust_ch
    path "virusdb.dmnd"
    publishDir "virusdb"

    """

    mmseqs createdb $params.query_file_aa DB
    mmseqs cluster --min-seq-id $params.mmseqs_minseqid --cov-mode 1 -c $params.mmseqs_cover DB DB_clu tmp
    mmseqs createsubdb DB_clu DB DB_clu_rep
    mmseqs convert2fasta DB_clu_rep DB_clu_rep.fasta
    diamond makedb --in DB_clu_rep.fasta -d virusdb
    rm -r tmp

    """

}

process hmmer {

    input:
    path x

    output:
    path "query_domains.hmmer"
    publishDir "query_domains"

    """

    hmmscan --cpu 2 --noali --notextw --qformat fasta --domtblout raw_domains.txt $params.phmms/*.hmm $x 1> /dev/null

    # Post-processing:
    # Merge overlapping query protein alignments, keep best (by bitscore)
    # For each best alignment, if e-value is under 0.1, report:
    # query overlap_start overlap_end best_start best_end best_bitscore best_i-Evalue model_acc model_name model_description

    grep -v "#" raw_domains.txt | \
    tr -s ' ' '\t' | \
    awk 'BEGIN{OFS="\t"}; {print \$4, \$18, \$19, \$13, \$14, \$2, \$1, \$23}' | \
    sort -k1,1 -k2,2n | \
    tee >(bedtools merge > tmp.out) | \
    bedtools cluster | \
    sort -k9,9n -k5,5nr | \
    sort -u -k9,9n | \
    bedtools intersect -a - -b tmp.out -wb -sorted | \
    awk 'BEGIN{OFS="\t"}; {if (\$4 < 0.1) print \$1, \$11, \$12, \$2, \$3, \$5, \$4, \$6, \$7, \$8}' \
    > query_domains.hmmer

    rm tmp.out

    """

}

process parse_ftp {

    output:
    path "*.ftp.txt"

    """

    while read line
        do
            assembly_ID=`echo \$line | sed 's/^.*\\///'`
            echo \$line > \${assembly_ID}.ftp.txt
        done < $params.ftp_file

    """

}

process download_assemblies {

    maxForks 8

    input:
    path x

    output:
    path "*genomic.fna.gz"

    """

    max_attempts=5

    count=0

    while read line
    do

    # Downloads and checks assembly file for corruption, re-attempts if md5 check fails
    md5check_function () {
        assembly_ID=`echo \$line | sed 's/^.*\\///'`
        wget -q \$line/\${assembly_ID}_genomic.fna.gz
        wget -q \$line/md5checksums.txt
        assemblyFile=`ls -1 *genomic.fna.gz`
        grep \$assemblyFile md5checksums.txt > \$assemblyFile.md5; rm md5checksums.txt
        status=`md5sum -c \$assemblyFile.md5 2>/dev/null | sed 's/.* //'`
        if [ "\$status" == FAILED ]
        then
                        if [ "\$count" == "\$max_attempts" ]
                        then
                                echo "\$assemblyFile FAILED md5check \$max_attempts times, exiting"; exit 1
                        else
                                echo "\$assemblyFile FAILED md5check"; rm \$assemblyFile*; count=\$count +1; md5check_function
                        fi
        else
                echo "\$assemblyFile PASSED md5check"; rm \$assemblyFile.md5
        fi
        }

    md5check_function

    done < $x

    """
}

process assembly_stats {

    input:
    path x

    output:
    path "assembly_stats.txt"

    """

    stats.sh in=$x format=3 addname= | grep -v n_scaffolds | sed 's/\\/.*\\///g; s/_genomic.fna.gz//' > assembly_stats.txt

    """

}

process diamond {

    maxForks 4

    input:
    path x

    output:
    path "*.dmnd.tsv", emit: tsv_ch
    path "*.gz*nsq"

    """

    db=$PWD/virusdb/virusdb.dmnd
    assembly=$x

    idle_fn () {
        if test -f "\$db"
            then
                return
            else
                sleep 0.5
                idle_fn
        fi
    }

    idle_fn

    diamond blastx \
    --$params.diamond_mode \
    --matrix $params.diamond_matrix \
    --masking seg \
    -d \$db \
    -q $x \
    -o matches.dmnd.tsv \
    -p $params.diamond_cpus \
    --outfmt 6 qseqid qstart qend qframe qlen sseqid sstart send slen evalue bitscore pident length mismatch gapopen &

    gunzip -c $x | makeblastdb -in - -out \${assembly/*\\//} -title \${assembly/*\\//} -dbtype nucl -parse_seqids &

    wait

    rm "\$(readlink -f $x)"

    """

}

process merge_and_extract {

    input:
    path x
    path y

    output:
    path "*.nonredundant.bed"
    path "*.fasta"

    """

    awk 'BEGIN{OFS="\t"}; {if(\$2<\$3) print \$1, \$2, \$3; else if(\$3<\$2) print \$1, \$3, \$2}' $x | \
    sort -k1,1 -k2,2n | bedtools merge -c 1 -o count > diamond-result.nonredundant.bed

    awk '{print \$1, \$2"-"\$3}' diamond-result.nonredundant.bed > batch.txt

    dbpath=\$(echo $y | cut -d ' ' -f1)
    dbpath=\$(readlink -f \$dbpath | sed 's/.nsq//g; s/\\.[0-9][0-9]\$//g')

    blastdbcmd -entry_batch batch.txt -db \$dbpath > result.fasta

    """

}

process intersect_domains {

    input:
    path x

    """

    domain_annotation=$PWD/query_domains/query_domains.hmmer

    idle_fn () {
        if test -f "\$domain_annotation"
            then
                return
            else
                sleep 0.5
                idle_fn
        fi
    }

    idle_fn

    awk 'BEGIN{OFS="\t"}; {if(\$2<\$3) print \$0; else if (\$3<\$2) print \$1, \$3, \$2, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13, \$14, \$15}' $x | \
    sort -k1,1 -k2,2n | \
    tee >(bedtools merge > tmp.out) | \
    bedtools cluster | \
    sort -k16,16n -k11,11nr | \
    sort -u -k16,16n | \
    bedtools intersect -a - -b tmp.out -loj -sorted | \
    awk 'BEGIN{OFS="\t"}; {print \$6, \$7, \$8, \$17, \$18, \$19, \$2, \$3, \$4, \$5, \$9, \$10, \$11, \$12, \$13, \$14, \$15}' | \
    sort -k1,1 -k2,2n | \
    bedtools intersect -a - -b \$domain_annotation -loj -sorted | \
    cut -f18 --complement > matches.dmnd.annot.tsv

    rm tmp.out

    """

}

// Workflow definition

workflow {

    // Build clustered DIAMOND query database from user supplied protein fasta
        build_db()

    // HMMER run on clustered queries
        hmmer (build_db.out.clust_ch)

    // Unpack user supplied ftp list and begin downloading assemblies
        fetched_assembly_files = parse_ftp() | flatten | download_assemblies

    // Sub-workflows to run on the downloaded assembly files
        assembly_stats(fetched_assembly_files)
        diamond(fetched_assembly_files) | merge_and_extract
        intersect_domains(diamond.out.tsv_ch)

}
