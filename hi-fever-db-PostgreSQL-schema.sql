-- hi-fever-db

-- PostgreSQL schema generation script

-- Dumped by pg_dump version 16.0
-- Maintainers: cormac.kinsella@evobio.eu, jose.ninobarreat@biology.ox.ac.uk

-- Instructions:
-- 1. Select a server, or create a new one (host = localhost, user = postgres)
-- 2. Create a new database within the server (e.g., hi-fever-db)
-- 3. Ensure sql data tables are located in a server accessible location
-- 4. Edit the import section of script to point to the HI-FEVER output files, lines 249-265)
-- 5. Right-click database name in the object explorer -> Query tool -> upload/paste this script and run it


-- MAIN SCRIPT --


-- Configure

SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', 'hifever_schema', false);
SET check_function_bodies = false;
SET xmloption = content;
SET client_min_messages = warning;
SET row_security = off;


-- Create schema

CREATE SCHEMA hifever_schema;

ALTER SCHEMA hifever_schema OWNER TO postgres;


-- Create tables

SET default_tablespace = '';

SET default_table_access_method = heap;

CREATE TABLE hifever_schema.assembly_statistics (
    scaffold_n bigint NOT NULL,
    contig_n bigint NOT NULL,
    scaffold_bp bigint NOT NULL,
    contig_bp bigint NOT NULL,
    gap_percent double precision NOT NULL,
    scaffold_n50 bigint NOT NULL,
    scaffold_l50 bigint NOT NULL,
    contig_n50 bigint NOT NULL,
    contig_l50 bigint NOT NULL,
    scaffold_n90 bigint NOT NULL,
    scaffold_l90 bigint NOT NULL,
    contig_n90 bigint NOT NULL,
    contig_l90 bigint NOT NULL,
    scaffold_max_length bigint NOT NULL,
    contig_max_length bigint NOT NULL,
    scaffold_n_greaterthan50k bigint NOT NULL,
    scaffold_percent_greaterthan50k double precision NOT NULL,
    gc_average double precision NOT NULL,
    gc_std double precision NOT NULL,
    assembly_id text NOT NULL
);

ALTER TABLE hifever_schema.assembly_statistics OWNER TO postgres;

CREATE TABLE hifever_schema.assembly_metadata (
    assembly_id text NOT NULL,
    assembly_accession text NOT NULL,
    bioproject text NOT NULL,
    biosample text NOT NULL,
    wgs_master text NOT NULL,
    refseq_category text NOT NULL,
    taxid text NOT NULL,
    species_taxid text NOT NULL,
    organism_name text NOT NULL,
    infraspecific_name text NOT NULL,
    isolate text NOT NULL,
    version_status text NOT NULL,
    assembly_level text NOT NULL,
    release_type text NOT NULL,
    genome_rep text NOT NULL,
    seq_rel_date text NOT NULL,
    asm_name text NOT NULL,
    asm_submitter text NOT NULL,
    gbrs_paired_asm text NOT NULL,
    paired_asm_comp text NOT NULL,
    ftp_path text NOT NULL,
    excluded_from_refseq text NOT NULL,
    relation_to_type_material text NOT NULL,
    asm_not_live_date text NOT NULL,
    assembly_type text NOT NULL,
    "group" text NOT NULL,
    genome_size text NOT NULL,
    genome_size_ungapped text NOT NULL,
    gc_percent text NOT NULL,
    replicon_count text NOT NULL,
    scaffold_count text NOT NULL,
    contig_count text NOT NULL,
    annotation_provider text NOT NULL,
    annotation_name text NOT NULL,
    annotation_date text NOT NULL,
    total_gene_count text NOT NULL,
    protein_coding_gene_count text NOT NULL,
    non_coding_gene_count text NOT NULL,
    pubmed_id text NOT NULL
);

ALTER TABLE hifever_schema.assembly_metadata OWNER TO postgres;

CREATE TABLE hifever_schema.best_forward_hits_pHMM (
    protein_id text NOT NULL,
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
    e_value_best_single double precision NOT NULL,
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
    domains_i_evalue_best_single double precision,
    model_acc text,
    model_name text,
    model_description text
);

ALTER TABLE hifever_schema.best_forward_hits_pHMM OWNER TO postgres;

CREATE TABLE hifever_schema.genewise (
    contig text NOT NULL,
    genomic_start bigint NOT NULL,
    genomic_end bigint NOT NULL,
    strand text NOT NULL,
    locus text NOT NULL,
    source_fasta text NOT NULL,
    bitscore double precision NOT NULL,
    query text NOT NULL,
    query_start bigint NOT NULL,
    query_end bigint NOT NULL,
    cdna_seq text NOT NULL,
    peptide_seq text NOT NULL,
    introns_predicted bigint NOT NULL,
    frameshifts_corrected bigint NOT NULL,
    inframe_stops_removed bigint NOT NULL
);

ALTER TABLE hifever_schema.genewise OWNER TO postgres;

CREATE TABLE hifever_schema.locus_to_assembly_junction_table (
    locus text NOT NULL,
    assembly_id text NOT NULL
);

ALTER TABLE hifever_schema.locus_to_assembly_junction_table OWNER TO postgres;

CREATE TABLE hifever_schema.predicted_orfs (
    locus text NOT NULL,
    coverage_of_locus_by_orf double precision NOT NULL,
    coverage_of_orf_by_locus double precision NOT NULL,
    orf_start bigint NOT NULL,
    orf_end bigint NOT NULL,
    strand text NOT NULL,
    orf_seq text NOT NULL
);

ALTER TABLE hifever_schema.predicted_orfs OWNER TO postgres;

CREATE TABLE hifever_schema.reciprocal_nr (
    query_locus text NOT NULL,
    subject_protein text NOT NULL,
    percent_identity double precision NOT NULL,
    "length" bigint NOT NULL,
    mismatches bigint NOT NULL,
    gapopens bigint NOT NULL,
    query_start bigint NOT NULL,
    query_end bigint NOT NULL,
    subject_start bigint NOT NULL,
    subject_end bigint NOT NULL,
    e_value double precision NOT NULL,
    bitscore double precision NOT NULL,
    subject_taxids text,
    subject_scientific_name text,
    subject_superkingdom text,
    subject_kingdom text,
    subject_phylum text,
    subject_title text NOT NULL
);

ALTER TABLE hifever_schema.reciprocal_nr OWNER TO postgres;

CREATE TABLE hifever_schema.reciprocal_rvdb (
    query_locus text NOT NULL,
    subject_protein text NOT NULL,
    percent_identity double precision NOT NULL,
    "length" bigint NOT NULL,
    mismatches bigint NOT NULL,
    gapopens bigint NOT NULL,
    query_start bigint NOT NULL,
    query_end bigint NOT NULL,
    subject_start bigint NOT NULL,
    subject_end bigint NOT NULL,
    e_value double precision NOT NULL,
    bitscore double precision NOT NULL,
    subject_taxids text,
    subject_scientific_name text,
    subject_superkingdom text,
    subject_kingdom text,
    subject_phylum text,
    subject_title text NOT NULL
);

ALTER TABLE hifever_schema.reciprocal_rvdb OWNER TO postgres;

CREATE TABLE hifever_schema.taxonomy (
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

ALTER TABLE hifever_schema.taxonomy OWNER TO postgres;


-- Import data

COPY assembly_statistics FROM '/HI-FEVER_OUTPUT_PATH/sql/assembly_stats.tsv' WITH DELIMITER E'\t' CSV;

COPY assembly_metadata FROM '/HI-FEVER_OUTPUT_PATH/sql/assembly_metadata.tsv' WITH DELIMITER E'\t' CSV;

COPY best_forward_hits_pHMM FROM '/HI-FEVER_OUTPUT_PATH/sql/matches.dmnd.annot.tsv' WITH DELIMITER E'\t' NULL '.' CSV;

COPY genewise FROM '/HI-FEVER_OUTPUT_PATH/sql/genewise.tsv' WITH DELIMITER E'\t' CSV;

COPY locus_to_assembly_junction_table FROM '/HI-FEVER_OUTPUT_PATH/sql/locus_assembly_map.tsv' WITH DELIMITER E'\t' CSV;

COPY predicted_orfs FROM '/HI-FEVER_OUTPUT_PATH/sql/predicted_ORFs.tsv' WITH DELIMITER E'\t' CSV;

COPY reciprocal_nr FROM '/HI-FEVER_OUTPUT_PATH/sql/reciprocal-nr-matches.dmnd.tsv' WITH DELIMITER E'\t' CSV;

COPY reciprocal_rvdb FROM '/HI-FEVER_OUTPUT_PATH/sql/reciprocal-rvdb-matches.dmnd.tsv' WITH DELIMITER E'\t' CSV;

COPY taxonomy FROM '/HI-FEVER_OUTPUT_PATH/sql/taxonomy_table.tsv' WITH DELIMITER E'\t' CSV;


-- Add primary key constraints

ALTER TABLE ONLY hifever_schema.assembly_statistics
    ADD CONSTRAINT assembly_statistics_pkey PRIMARY KEY (assembly_id);

ALTER TABLE ONLY hifever_schema.genewise
    ADD CONSTRAINT genewise_pkey PRIMARY KEY (locus);

ALTER TABLE ONLY hifever_schema.locus_to_assembly_junction_table
    ADD CONSTRAINT locus_to_assembly_junction_table_pkey PRIMARY KEY (locus);

ALTER TABLE ONLY hifever_schema.taxonomy
    ADD CONSTRAINT species_taxid_pkey PRIMARY KEY (species_taxid);


-- Add foreign key constraints

ALTER TABLE ONLY hifever_schema.best_forward_hits_pHMM
    ADD CONSTRAINT best_forward_hits_pHMM_locus_fkey FOREIGN KEY (locus) REFERENCES hifever_schema.locus_to_assembly_junction_table(locus) NOT VALID;

ALTER TABLE ONLY hifever_schema.locus_to_assembly_junction_table
    ADD CONSTRAINT link_junction_table_to_assembly_stats_fkey FOREIGN KEY (assembly_id) REFERENCES hifever_schema.assembly_statistics(assembly_id) NOT VALID;

ALTER TABLE ONLY hifever_schema.assembly_metadata
    ADD CONSTRAINT link_metadata_to_assembly_stats_fkey FOREIGN KEY (assembly_id) REFERENCES hifever_schema.assembly_statistics(assembly_id) NOT VALID;

ALTER TABLE ONLY hifever_schema.assembly_metadata
    ADD CONSTRAINT link_metadata_to_full_taxonomy_fkey FOREIGN KEY (species_taxid) REFERENCES hifever_schema.taxonomy(species_taxid) NOT VALID;

ALTER TABLE ONLY hifever_schema.locus_to_assembly_junction_table
    ADD CONSTRAINT link_to_genewise_fkey FOREIGN KEY (locus) REFERENCES hifever_schema.genewise(locus) NOT VALID;

ALTER TABLE ONLY hifever_schema.predicted_orfs
    ADD CONSTRAINT predicted_orfs_locus_fkey FOREIGN KEY (locus) REFERENCES hifever_schema.locus_to_assembly_junction_table(locus) NOT VALID;

ALTER TABLE ONLY hifever_schema.reciprocal_nr
    ADD CONSTRAINT reciprocal_nr_query_locus_fkey FOREIGN KEY (query_locus) REFERENCES hifever_schema.locus_to_assembly_junction_table(locus) NOT VALID;

ALTER TABLE ONLY hifever_schema.reciprocal_nr
    ADD CONSTRAINT link_recip_nr_to_full_taxonomy_fkey FOREIGN KEY (subject_taxids) REFERENCES hifever_schema.taxonomy(species_taxid) NOT VALID;

ALTER TABLE ONLY hifever_schema.reciprocal_rvdb
    ADD CONSTRAINT reciprocal_rvdb_query_locus_fkey FOREIGN KEY (query_locus) REFERENCES hifever_schema.locus_to_assembly_junction_table(locus) NOT VALID;

ALTER TABLE ONLY hifever_schema.reciprocal_rvdb
    ADD CONSTRAINT link_recip_rvdb_to_full_taxonomy_fkey FOREIGN KEY (subject_taxids) REFERENCES hifever_schema.taxonomy(species_taxid) NOT VALID;
