#!/home/user/nextflow

// Syntax version

nextflow.enable.dsl=2

// Script parameters

params.ftp_file = "$PWD/ftp_list.txt"
params.query_file_aa = "$PWD/protein_query.fasta"
params.phmms = "$PWD/profiles.hmm"
params.mmseqs_minseqid = "0.95"
params.mmseqs_cover = "0.90"
params.diamond_mode = "very-sensitive"
params.diamond_matrix = "BLOSUM62"
params.diamond_cpus = "12"

// Define workflow processes

process build_db {

    input:
    path x

    output:
    path "DB_clu_rep.fasta", emit: clust_ch
    path "virusdb.dmnd"
    publishDir "virusdb"

    """

    mmseqs createdb $x DB
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
    path y

    output:
    path "query_domains.hmmer"

    """

    hmmscan --noali --notextw --qformat fasta --domtblout raw_out.txt $x $y
    # --cpu 10. Or parameter?

    # POST PROCESSING
    # raw_out.txt > query_domains.hmmer

    """

}

process parse_ftp {

    input:
    path x

    output:
    path "*.ftp.txt"

    """

    while read line
    do

    assembly_ID=`echo \$line | sed 's/^.*\\///'`

    echo \$line > \${assembly_ID}.ftp.txt

    done < $x

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
    path "*.dmnd.tsv"
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

process bedtools_extract {

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

// Workflow definition

workflow {

    // Build clustered diamond query database from user supplied protein fasta
        def db_ch = Channel.fromPath(params.query_file_aa)
        build_db(db_ch)
        def profiles_ch = Channel.fromPath(params.phmms)
        hmmer (profiles_ch, build_db.out.clust_ch)

    // Unpack user supplied ftp list and begin downloading assemblies
        def ftp_ch = Channel.fromPath(params.ftp_file)
        fetched_assembly_files = parse_ftp(ftp_ch).flatten() | download_assemblies

    // Independent sub-workflows to run on the downloaded assembly files
        assembly_stats(fetched_assembly_files)
        diamond(fetched_assembly_files) | bedtools_extract

}
