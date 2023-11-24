-- hi-fever-db

-- PostgreSQL schema generation script

-- Dumped by pg_dump version 16.0
-- Maintainers: cormac.kinsella@evobio.eu, jose.ninobarreat@biology.ox.ac.uk

-- Instructions:
-- 1. Select a server, or create a new one (host = localhost, user = postgres)
-- 2. Create a new database within the server (e.g., hi-fever-db)
-- 3. Ensure sql data tables are located in a server accessible location (& edit import section of script accordingly, lines 246-264)
-- 4. Right-click database name in the object explorer -> Query tool -> upload/paste this script and run it


-- MAIN SCRIPT --


-- Configure

SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', 'hi-fever-schema', false);
SET check_function_bodies = false;
SET xmloption = content;
SET client_min_messages = warning;
SET row_security = off;


-- Create schema

CREATE SCHEMA "hi-fever-schema";

ALTER SCHEMA "hi-fever-schema" OWNER TO postgres;


-- Create tables

SET default_tablespace = '';

SET default_table_access_method = heap;

CREATE TABLE "hi-fever-schema"."assembly-statistics" (
    "scaffold_N" bigint NOT NULL,
    "contig_N" bigint NOT NULL,
    scaffold_bp bigint NOT NULL,
    contig_bp bigint NOT NULL,
    gap_percent double precision NOT NULL,
    "scaffold_N50" bigint NOT NULL,
    "scaffold_L50" bigint NOT NULL,
    "contig_N50" bigint NOT NULL,
    "contig_L50" bigint NOT NULL,
    "scaffold_N90" bigint NOT NULL,
    "scaffold_L90" bigint NOT NULL,
    "contig_N90" bigint NOT NULL,
    "contig_L90" bigint NOT NULL,
    scaffold_max_length bigint NOT NULL,
    contig_max_length bigint NOT NULL,
    "scaffold_N_greaterthan50k" bigint NOT NULL,
    "scaffold_percent_greaterthan50K" double precision NOT NULL,
    gc_average double precision NOT NULL,
    gc_std double precision NOT NULL,
    "assembly-ID" text NOT NULL
);

ALTER TABLE "hi-fever-schema"."assembly-statistics" OWNER TO postgres;

CREATE TABLE "hi-fever-schema"."assembly-metadata" (
    "assembly-ID" text NOT NULL,
    "assembly_accession" text NOT NULL,
    "bioproject" text NOT NULL,
    "biosample" text NOT NULL,
    "wgs_master" text NOT NULL,
    "refseq_category" text NOT NULL,
    "taxid" text NOT NULL,
    "species_taxid" text NOT NULL,
    "organism_name" text NOT NULL,
    "infraspecific_name" text NOT NULL,
    "isolate" text NOT NULL,
    "version_status" text NOT NULL,
    "assembly_level" text NOT NULL,
    "release_type" text NOT NULL,
    "genome_rep" text NOT NULL,
    "seq_rel_date" text NOT NULL,
    "asm_name" text NOT NULL,
    "asm_submitter" text NOT NULL,
    "gbrs_paired_asm" text NOT NULL,
    "paired_asm_comp" text NOT NULL,
    "ftp_path" text NOT NULL,
    "excluded_from_refseq" text NOT NULL,
    "relation_to_type_material" text NOT NULL,
    "asm_not_live_date" text NOT NULL,
    "assembly_type" text NOT NULL,
    "group" text NOT NULL,
    "genome_size" text NOT NULL,
    "genome_size_ungapped" text NOT NULL,
    "gc_percent" text NOT NULL,
    "replicon_count" text NOT NULL,
    "scaffold_count" text NOT NULL,
    "contig_count" text NOT NULL,
    "annotation_provider" text NOT NULL,
    "annotation_name" text NOT NULL,
    "annotation_date" text NOT NULL,
    "total_gene_count" text NOT NULL,
    "protein_coding_gene_count" text NOT NULL,
    "non_coding_gene_count" text NOT NULL,
    "pubmed_id" text NOT NULL
);

ALTER TABLE "hi-fever-schema"."assembly-metadata" OWNER TO postgres;

CREATE TABLE "hi-fever-schema"."best-forward-hits-pHMM" (
    "protein_ID" text NOT NULL,
    protein_start bigint NOT NULL,
    protein_end bigint NOT NULL,
    locus text NOT NULL,
    contig text NOT NULL,
    contig_start_merged_aligns bigint NOT NULL,
    contig_end_merged_aligns bigint NOT NULL,
    contig_start_best_single bigint NOT NULL,
    contig_end_best_single bigint NOT NULL,
    frame_best_single text NOT NULL,
    contig_length_nt bigint NOT NULL,
    protein_length_aa_best_single bigint NOT NULL,
    "e-value_best_single" double precision NOT NULL,
    bitscore_best_single double precision NOT NULL,
    percent_identity_best_single double precision NOT NULL,
    align_length_aa_best_single bigint NOT NULL,
    mismatches_best_single bigint NOT NULL,
    gapopen_best_single bigint NOT NULL,
    domains_merged_start bigint NOT NULL,
    domains_merged_end bigint NOT NULL,
    domains_start_best_single bigint,
    domains_end_best_single bigint,
    domains_bitscore_best_single double precision,
    "domains_i-evalue_best_single" double precision,
    model_acc text,
    model_name text,
    model_description text
);

ALTER TABLE "hi-fever-schema"."best-forward-hits-pHMM" OWNER TO postgres;

CREATE TABLE "hi-fever-schema".genewise (
    contig text NOT NULL,
    genomic_start bigint NOT NULL,
    genomic_end bigint NOT NULL,
    strand text NOT NULL,
    locus text NOT NULL,
    "sourceFASTA" text NOT NULL,
    bitscore double precision NOT NULL,
    query text NOT NULL,
    query_start bigint NOT NULL,
    query_end bigint NOT NULL,
    cdna_seq text NOT NULL,
    peptide_seq text NOT NULL,
    introns_predicted bigint NOT NULL,
    frameshifts_corrected bigint NOT NULL,
    "inframe_STOPs_removed" bigint NOT NULL
);

ALTER TABLE "hi-fever-schema".genewise OWNER TO postgres;

CREATE TABLE "hi-fever-schema"."locus-to-assembly-junction-table" (
    locus text NOT NULL,
    "assembly-ID" text NOT NULL
);

ALTER TABLE "hi-fever-schema"."locus-to-assembly-junction-table" OWNER TO postgres;

CREATE TABLE "hi-fever-schema"."predicted-orfs" (
    locus text NOT NULL,
    "coverage_of_locus_by_ORF" double precision NOT NULL,
    "coverage_of_ORF_by_locus" double precision NOT NULL,
    "ORF_start" bigint NOT NULL,
    "ORF_end" bigint NOT NULL,
    strand text NOT NULL,
    "ORF_seq" text NOT NULL
);

ALTER TABLE "hi-fever-schema"."predicted-orfs" OWNER TO postgres;

CREATE TABLE "hi-fever-schema"."reciprocal-nr" (
    query_locus text NOT NULL,
    subject_protein text NOT NULL,
    percent_identity double precision NOT NULL,
    length bigint NOT NULL,
    mismatches bigint NOT NULL,
    gapopens bigint NOT NULL,
    query_start bigint NOT NULL,
    query_end bigint NOT NULL,
    subject_start bigint NOT NULL,
    subject_end bigint NOT NULL,
    "e-value" double precision NOT NULL,
    bitscore double precision NOT NULL,
    subject_taxids text,
    subject_scientific_name text,
    subject_superkingdom text,
    subject_kingdom text,
    subject_phylum text,
    subject_title text NOT NULL
);

ALTER TABLE "hi-fever-schema"."reciprocal-nr" OWNER TO postgres;

CREATE TABLE "hi-fever-schema"."reciprocal-rvdb" (
    query_locus text NOT NULL,
    subject_protein text NOT NULL,
    percent_identity double precision NOT NULL,
    length bigint NOT NULL,
    mismatches bigint NOT NULL,
    gapopens bigint NOT NULL,
    query_start bigint NOT NULL,
    query_end bigint NOT NULL,
    subject_start bigint NOT NULL,
    subject_end bigint NOT NULL,
    "e-value" double precision NOT NULL,
    bitscore double precision NOT NULL,
    subject_taxids text,
    subject_scientific_name text,
    subject_superkingdom text,
    subject_kingdom text,
    subject_phylum text,
    subject_title text NOT NULL
);

ALTER TABLE "hi-fever-schema"."reciprocal-rvdb" OWNER TO postgres;

CREATE TABLE "hi-fever-schema"."taxonomy" (
    species_taxid text NOT NULL,
    superkingdom text NOT NULL,
    kingdom text NOT NULL,
    phylum text NOT NULL,
    class text NOT NULL,
    "order" text NOT NULL,
    family text NOT NULL,
    genus text NOT NULL,
    species text NOT NULL
);

ALTER TABLE "hi-fever-schema"."taxonomy" OWNER TO postgres;


-- Import data

COPY "assembly-statistics" FROM 'C:/Program Files/PostgreSQL/16/data/sql/assembly_stats.tsv' WITH DELIMITER E'\t' CSV;

COPY "assembly-metadata" FROM 'C:/Program Files/PostgreSQL/16/data/sql/assembly_metadata.tsv' WITH DELIMITER E'\t' CSV;

COPY "best-forward-hits-pHMM" FROM 'C:/Program Files/PostgreSQL/16/data/sql/matches.dmnd.annot.tsv' WITH DELIMITER E'\t' NULL '.' CSV;

COPY "genewise" FROM 'C:/Program Files/PostgreSQL/16/data/sql/genewise.tsv' WITH DELIMITER E'\t' CSV;

COPY "locus-to-assembly-junction-table" FROM 'C:/Program Files/PostgreSQL/16/data/sql/locus_assembly_map.tsv' WITH DELIMITER E'\t' CSV;

COPY "predicted-orfs" FROM 'C:/Program Files/PostgreSQL/16/data/sql/predicted_ORFs.tsv' WITH DELIMITER E'\t' CSV;

COPY "reciprocal-nr" FROM 'C:/Program Files/PostgreSQL/16/data/sql/reciprocal-nr-matches.dmnd.tsv' WITH DELIMITER E'\t' CSV;

COPY "reciprocal-rvdb" FROM 'C:/Program Files/PostgreSQL/16/data/sql/reciprocal-rvdb-matches.dmnd.tsv' WITH DELIMITER E'\t' CSV;

COPY "taxonomy" FROM 'C:/Program Files/PostgreSQL/16/data/sql/taxonomy_table.tsv' WITH DELIMITER E'\t' CSV;


-- Add primary key constraints

ALTER TABLE ONLY "hi-fever-schema"."assembly-statistics"
    ADD CONSTRAINT "assembly-statistics_pkey" PRIMARY KEY ("assembly-ID");

ALTER TABLE ONLY "hi-fever-schema".genewise
    ADD CONSTRAINT genewise_pkey PRIMARY KEY (locus);

ALTER TABLE ONLY "hi-fever-schema"."locus-to-assembly-junction-table"
    ADD CONSTRAINT "locus-to-assembly-junction-table_pkey" PRIMARY KEY (locus);

ALTER TABLE ONLY "hi-fever-schema"."taxonomy"
    ADD CONSTRAINT "species-taxid_pkey" PRIMARY KEY (species_taxid);


-- Add foreign key constraints

ALTER TABLE ONLY "hi-fever-schema"."best-forward-hits-pHMM"
    ADD CONSTRAINT "best-forward-hits-pHMM_locus_fkey" FOREIGN KEY (locus) REFERENCES "hi-fever-schema"."locus-to-assembly-junction-table"(locus) NOT VALID;

ALTER TABLE ONLY "hi-fever-schema"."locus-to-assembly-junction-table"
    ADD CONSTRAINT "link-junction-table-to-assembly-stats_fkey" FOREIGN KEY ("assembly-ID") REFERENCES "hi-fever-schema"."assembly-statistics"("assembly-ID") NOT VALID;

ALTER TABLE ONLY "hi-fever-schema"."assembly-metadata"
    ADD CONSTRAINT "link-metadata-to-assembly-stats_fkey" FOREIGN KEY ("assembly-ID") REFERENCES "hi-fever-schema"."assembly-statistics"("assembly-ID") NOT VALID;

ALTER TABLE ONLY "hi-fever-schema"."assembly-metadata"
    ADD CONSTRAINT "link-metadata-to-full-taxonomy_fkey" FOREIGN KEY (species_taxid) REFERENCES "hi-fever-schema"."taxonomy"(species_taxid) NOT VALID;

ALTER TABLE ONLY "hi-fever-schema"."locus-to-assembly-junction-table"
    ADD CONSTRAINT "link-to-genewise_fkey" FOREIGN KEY (locus) REFERENCES "hi-fever-schema".genewise(locus) NOT VALID;

ALTER TABLE ONLY "hi-fever-schema"."predicted-orfs"
    ADD CONSTRAINT "predicted-orfs_locus_fkey" FOREIGN KEY (locus) REFERENCES "hi-fever-schema"."locus-to-assembly-junction-table"(locus) NOT VALID;

ALTER TABLE ONLY "hi-fever-schema"."reciprocal-nr"
    ADD CONSTRAINT "reciprocal-nr_query_locus_fkey" FOREIGN KEY (query_locus) REFERENCES "hi-fever-schema"."locus-to-assembly-junction-table"(locus) NOT VALID;

ALTER TABLE ONLY "hi-fever-schema"."reciprocal-nr"
    ADD CONSTRAINT "link-recip-nr-to-full-taxonomy_fkey" FOREIGN KEY (subject_taxids) REFERENCES "hi-fever-schema"."taxonomy"(species_taxid) NOT VALID;

ALTER TABLE ONLY "hi-fever-schema"."reciprocal-rvdb"
    ADD CONSTRAINT "reciprocal-rvdb_query_locus_fkey" FOREIGN KEY (query_locus) REFERENCES "hi-fever-schema"."locus-to-assembly-junction-table"(locus) NOT VALID;

ALTER TABLE ONLY "hi-fever-schema"."reciprocal-rvdb"
    ADD CONSTRAINT "link-recip-rvdb-to-full-taxonomy_fkey" FOREIGN KEY (subject_taxids) REFERENCES "hi-fever-schema"."taxonomy"(species_taxid) NOT VALID;
