#!/home/user/nextflow

// Syntax version

nextflow.enable.dsl=2

// Script parameters

params.ftp_file = "ftp_list.txt"
params.query_file_aa = "protein_query.fasta"
params.phmms = "domains"
params.mmseqs_minseqid = "0.95"
params.mmseqs_cover = "0.90"
params.diamond_mode = "very-sensitive"
params.diamond_matrix = "BLOSUM62"
params.diamond_cpus = "12"
params.interval = "1000"
params.flank = "3000"
params.reciprocal_db = "nr_clustered.dmnd"
params.orf_size_nt = "150"
params.genewise_matrix = "BLOSUM62"
params.stop_task = "remove"

// Define workflow processes

process build_db {

    input:
    path queries

    output:
    path "DB_clu_rep.fasta", emit: clust_ch
    path "virusdb.dmnd"
    publishDir "virusdb"

    """

    mmseqs createdb $queries DB
    mmseqs cluster --min-seq-id $params.mmseqs_minseqid --cov-mode 1 -c $params.mmseqs_cover DB DB_clu tmp
    mmseqs createsubdb DB_clu DB DB_clu_rep
    mmseqs convert2fasta DB_clu_rep DB_clu_rep.fasta
    diamond makedb --in DB_clu_rep.fasta -d virusdb
    rm -r tmp

    """

}

process hmmer {

    input:
    path profile_dir
    path clustered_fasta

    output:
    path "query_domains.hmmer"
    publishDir "query_domains"

    """

    hmmscan --cpu 4 --noali --notextw --qformat fasta --domtblout raw_domains.txt $profile_dir/*.hmm $clustered_fasta 1> /dev/null

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
    awk 'BEGIN{OFS="\t"}; {if (\$4 < 0.1) print \$1, \$11, \$12, \$2, \$3, \$5, \$4, \$6, \$7, \$8}' > \
    query_domains.hmmer

    rm tmp.out

    """

}

process parse_ftp {

    input:
    path ftp_input

    output:
    path "*.ftp.txt"

    """

    while read line
        do
            assemblyID=`echo \$line | sed 's/^.*\\///'`
            echo \$line > \${assemblyID}.ftp.txt
        done < $ftp_input

    """

}

process download_assemblies {

    maxForks 8

    input:
    path ftp_dir

    output:
    path "*genomic.fna.gz"

    """

    max_attempts=5

    count=0

    while read line
    do

    # Downloads and checks assembly file for corruption, re-attempts if md5 check fails
    md5check_function () {
        assemblyFile="\$(echo "\$line" | sed 's/^.*\\///')_genomic.fna.gz"
        wget -q "\$line/\$assemblyFile" "\$line/md5checksums.txt"
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

    done < $ftp_dir

    """
}

process assembly_stats {

    input:
    path assembly

    output:
    path "assembly_stats.txt"

    """

    stats.sh in=$assembly format=3 addname= | \
    grep -v n_scaffolds | \
    sed 's/\\/.*\\///g; s/_genomic.fna.gz//' > \
    assembly_stats.txt

    """

}

process diamond {

    maxForks 4

    input:
    path assembly

    output:
    path "*.dmnd.tsv"
    path "*.gz*nsq"

    """

    db=$PWD/virusdb/virusdb.dmnd
    assembly=$assembly

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
    -q \$assembly \
    -o matches.dmnd.tsv \
    -p $params.diamond_cpus \
    --outfmt 6 qseqid qstart qend qframe qlen sseqid sstart send slen evalue bitscore pident length mismatch gapopen &

    gunzip -c \$assembly | makeblastdb -in - -out \${assembly/*\\//} -title \${assembly/*\\//} -dbtype nucl -parse_seqids &

    wait

    rm "\$(readlink -f \$assembly)"

    """

}

process intersect_domains_merge_extract {

    input:
    path diamond_tsv
    path assembly_nsq_db

    output:
    path "*_strict.fasta", emit: strict_ch
    path "*_context.fasta", emit: context_ch
    path "diamond-result.nonredundant.bed", emit: nr_bed_ch
    path "matches.dmnd.annot.tsv", emit: annot_tsv_ch

    """

    # Wait for domain annotation file

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

    # Intersect domains & produce nonredundant BED:
    # Converts DIAMOND tsv to ascending assembly coordinate ranges, sorts to BED compatibility (contig and start position).
    # Calculates maximal strictly overlapping coordinate ranges and stores in temp file.
    # Adds a column to main tsv denoting which coordinates are in mergable overlapping clusters
    # Sorts best bitscore alignment per cluster to the top and removes others
    # Intersects this with the maximal range temp file
    # Converts to a subject oriented bed (i.e., protein hit and coordinates).
    # Intersects with domain coordinate annotation file, reports:
    # sseqid_(protein) sstart send qseqid_best qstart_overlap qend_overlap qstart_best qend_best qframe_best qlen_best slen \
    # eval bitscore pident len mismatch gapopen domain_overlap_start domain_overlap_end best_start best_end best_bitscore best_i-Evalue \
    # model_acc model_name model_description

    awk 'BEGIN{OFS="\t"}; {if(\$2<\$3) print \$0; else if (\$3<\$2) print \$1, \$3, \$2, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13, \$14, \$15}' $diamond_tsv | \
    sort -k1,1 -k2,2n | \
    tee >(bedtools merge > diamond-result.nonredundant.bed) | \
    bedtools cluster | \
    sort -k16,16n -k11,11nr | \
    sort -u -k16,16n | \
    bedtools intersect -a - -b diamond-result.nonredundant.bed -loj -sorted | \
    awk 'BEGIN{OFS="\t"}; {print \$6, \$7, \$8, \$17, \$18, \$19, \$2, \$3, \$4, \$5, \$9, \$10, \$11, \$12, \$13, \$14, \$15}' | \
    sort -k1,1 -k2,2n | \
    bedtools intersect -a - -b \$domain_annotation -loj -sorted | \
    cut -f18 --complement > \
    matches.dmnd.annot.tsv

    # Prepare variables

    dbpath=\$(readlink -f \$(echo $assembly_nsq_db | cut -d ' ' -f1) | sed 's/.nsq//g; s/\\.[0-9][0-9]\$//g')
    filename=\$(echo \$dbpath | sed 's/\\.gz//g; s/\\/.*\\///g')

    # First coordinate range extraction (strictly overlapping alignments)

    awk '{print \$1, \$2"-"\$3}' diamond-result.nonredundant.bed | \
    blastdbcmd -entry_batch - -db \$dbpath > \
    "\${filename}_strict.fasta"

    # Second coordinate range extraction (allow interval and add flanks)

    bedtools merge -d $params.interval -i diamond-result.nonredundant.bed | \
    awk -v flank=$params.flank '{if(\$2-flank < 1) print \$1, 1"-"\$3+flank; else print \$1, \$2-flank"-"\$3+flank}' | \
    blastdbcmd -entry_batch - -db \$dbpath > \
    "\${filename}_context.fasta"

    rm \$dbpath*

    """

}

process orf_extract {

    input:
    path context_fasta
    path nr_bed

    """

    # Checks input file for content, if yes:
    # Extracts ORFs between START and STOP codons. STOP is NOT included in reported ORF coordinates.
    # Converts ORF output to the original genomic loci in ascending BED format, with predicted protein sequence.

    if [ -s $context_fasta ];
        then
            cat $context_fasta | \
            sed 's/:/-/' | \
            getorf -sequence /dev/stdin -outseq /dev/stdout -minsize $params.orf_size_nt -find 1 2>/dev/null | \
            seqtk seq - | \
            sed -r 's/] .*//; s/_[0-9]+ \\[/ /; s/>//; s/-/ /g' | \
            tr -s ' ' '\t' | \
            paste -sd '\t\n' | \
            awk 'BEGIN{OFS="\t"}; {if (\$4 < \$5) print \$1, \$2+\$4-1, \$2+\$5-1, "+", \$6; else print \$1, \$2+\$5-1, \$2+\$4-1, "-", \$6}' | \
            sort -k1,1 -k2,2n > \
            context_ORFs.txt
        else
            touch context_ORFs.txt
    fi

    # Intersect ORFs with strictly overlapping features, reports:
    # contig strict_feature_start strict_feature_end coverage_of_feature_by_ORF coverage_of_ORF_by_feature ORF_start ORF_end ORF_strand ORF_seq

    bedtools intersect -a $nr_bed -b context_ORFs.txt -wo -sorted | \
    awk '{print \$1, \$2, \$3, \$9/(\$3-\$2), (\$3-\$2)/(\$6-\$5), \$5, \$6, \$7, \$8}' > \
    intersected_ORFs.txt

    """

}

process genewise {

    input:
    path annotated_tsv
    path strict_fasta

    """

    # Check for file content

    if [ -s $strict_fasta ];

        then

            protein_db=$PWD/virusdb/DB_clu_rep.fasta
            export WISECONFIGDIR="\$CONDA_PREFIX/share/wise2/wisecfg"
            mkdir wise_tmp

            # Reformat and extract nucleotide FASTA, one per file

            sed 's/ .*//' < $strict_fasta > wise_tmp/temp.fa
            grep ">" wise_tmp/temp.fa | sed 's/>//' > wise_tmp/nuc_headers

            while read line
                do
                    echo \$line | \
                    seqtk subseq wise_tmp/temp.fa - > \
                    wise_tmp/\$line
                done < wise_tmp/nuc_headers

            # Generate query-target pairing file & protein accessions to extract

            cut -f1,4-6 $annotated_tsv | \
            uniq | \
            tee >(cut -f1 | sort | uniq > wise_tmp/prot_headers) | \
            awk 'BEGIN{OFS="\t"} {print \$1,\$2":"\$3"-"\$4}' > \
            wise_tmp/matched_pairs

            # Extract protein FASTA, one per file

            while read line
                do
                    echo \$line | \
                    seqtk subseq \$protein_db - > \
                    wise_tmp/\$line
                done < wise_tmp/prot_headers

            # GeneWise operations

            while read line
                do
                    query=\$(echo \$line | cut -f1 -d ' ')
                    target=\$(echo \$line | cut -f2 -d ' ')
                    genewise wise_tmp/\$query wise_tmp/\$target -both -matrix "$params.genewise_matrix".bla -sum -pep -cdna -divide DIVIDE_STRING -silent | \
                    grep -v ">\\|Bits   Query" | \
                    awk '/^[-0-9]/ {printf("%s%s\\t",(N>0?"\\n":""),\$0);N++;next;} {printf("%s",\$0);} END {printf("\\n");}' | \
                    sed 's/DIVIDE_STRING/\t/g' | \
                    tr -s '\t' | \
                    tr -s ' ' '\t' | \
                    sort -k5,5 -k1,1nr | \
                    sort -u -k5,5 >> \
                    genewise
                done < wise_tmp/matched_pairs

            # Cleanup

            rm -r wise_tmp

            # Post-processing

            python $PWD/scripts/stop_convert_and_count.py --task $params.stop_task --file genewise > genewise_processed.txt
            rm "\$(readlink -f genewise)"

        else

            touch genewise_processed.txt

    fi

    """

}

process reciprocal_diamond {

    input:
    path strict_fastas
    path reciprocal_db

    output:
    path "*.dmnd.tsv"

    """

    cat $strict_fastas > all_seqs.fasta

    diamond blastx \
    --$params.diamond_mode \
    --matrix $params.diamond_matrix \
    --masking seg \
    -d $reciprocal_db \
    -q all_seqs.fasta \
    -o reciprocal-matches.dmnd.tsv \
    -e 1e-5 \
    -k 20 \
    --outfmt 6 qseqid qstart qend qframe qlen sseqid sstart send slen evalue bitscore pident length mismatch gapopen

    """

}

// Workflow definition

workflow {

    // Build clustered DIAMOND query database from user supplied protein fasta
        def query_ch = Channel.fromPath(params.query_file_aa)
        build_db(query_ch)

    // HMMER run on clustered queries
        def profiles_ch = Channel.fromPath(params.phmms, type: 'dir')
        hmmer (profiles_ch, build_db.out.clust_ch)

    // Unpack user supplied ftp list and begin downloading assemblies
        def ftp_ch = Channel.fromPath(params.ftp_file)
        fetched_assembly_files = parse_ftp(ftp_ch) | flatten | download_assemblies

    // Analyse downloaded assembly files
        assembly_stats(fetched_assembly_files)
        diamond(fetched_assembly_files) | intersect_domains_merge_extract

    // Extract ORFs that overlap DIAMOND hits (extending into flanks)
        orf_extract (intersect_domains_merge_extract.out.context_ch, intersect_domains_merge_extract.out.nr_bed_ch)

    // Frameshift and STOP aware reconstruction of EVE sequences
        genewise (intersect_domains_merge_extract.out.annot_tsv_ch, intersect_domains_merge_extract.out.strict_ch)

    // Reciprocal DIAMOND
        def reciprocal_db_ch = Channel.fromPath(params.reciprocal_db)
        collected_strict_fastas = intersect_domains_merge_extract.out.strict_ch.collect()
        reciprocal_diamond (collected_strict_fastas, reciprocal_db_ch)

}
