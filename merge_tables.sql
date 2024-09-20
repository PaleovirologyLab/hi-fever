--Wenyang Lyu, Jose Gabriel Nino Barreat, Emma Harding
--Date: 20/Sep/2024

-- Set schema to use

SET search_path TO test_schema;

--Add database names to reciprocal tables

ALTER TABLE reciprocal_rvdb
ADD database TEXT DEFAULT 'RVDB';

ALTER TABLE reciprocal_nr
ADD database TEXT DEFAULT 'nr';

--Merge all reciprocal matches and get subject sequence taxonomic information

CREATE TABLE all_reciprocal_matches AS
SELECT query_locus, "database", subject_protein,
ARRAY_TO_STRING(REGEXP_MATCHES(subject_title, '.+?\|.+?\|.+?\|.+?\|.+?\|(.+?).\[.+?\]$'), ', ') AS subject_title,
percent_identity, "length", mismatches, gapopens, e_value, bitscore, subject_taxids, superkingdom, kingdom,
phylum, "class", "order", "family", genus, species
FROM reciprocal_rvdb
LEFT JOIN taxonomy
ON taxonomy.species_taxid = reciprocal_rvdb.subject_taxids
UNION ALL
SELECT query_locus, "database", subject_protein,
ARRAY_TO_STRING(REGEXP_MATCHES(subject_title, '.+?\s(.+?)\['), ', ') AS subject_title,
percent_identity, "length", mismatches, gapopens, e_value, bitscore, subject_taxids, superkingdom, kingdom,
phylum, "class", "order", "family", genus, species
FROM reciprocal_nr
LEFT JOIN taxonomy
ON taxonomy.species_taxid = reciprocal_nr.subject_taxids;

--Create the large merged table

CREATE TABLE merged AS (
SELECT
    all_reciprocal_matches.query_locus as EVE_locus,
    all_reciprocal_matches.percent_identity as blast_identity,
    all_reciprocal_matches.e_value as blast_evalue,
    all_reciprocal_matches.bitscore as blast_bitscore,
    all_reciprocal_matches.subject_taxids as hit_taxid,
    all_reciprocal_matches.subject_title as hit_title,
    all_reciprocal_matches."database",
    all_reciprocal_matches."length" as hit_length,
    all_reciprocal_matches.kingdom as hit_kingdom,
    all_reciprocal_matches.phylum as hit_phylum,
    all_reciprocal_matches.class as hit_class,
    all_reciprocal_matches.order as hit_order,
    all_reciprocal_matches.family as hit_family,
    all_reciprocal_matches.genus as hit_genus,
    assembly_metadata.assembly_id,
    assembly_metadata.genome_size_ungapped as host_genome_size,
    assembly_metadata.organism_name as host_species,
    assembly_metadata.taxid as host_taxid,
    assembly_metadata.biosample as host_biosample,
    assembly_metadata."group" as host_group,
    genewise.strand,
    genewise.introns_predicted,
    genewise.frameshifts_corrected,
    genewise.inframe_stops_removed,
    genewise.query as genewise_model,
    genewise.peptide_seq as reconstructed_peptide,
    genewise.bitscore as genewise_bitscore,
    taxonomy."class" as host_class,
    taxonomy."order" as host_order,
    taxonomy.family as host_family,
    taxonomy.genus as host_genus
FROM
    all_reciprocal_matches
    JOIN locus_to_assembly_junction_table
    ON all_reciprocal_matches.query_locus = locus_to_assembly_junction_table.locus
    JOIN assembly_metadata
    ON assembly_metadata.assembly_id = locus_to_assembly_junction_table.assembly_id
    JOIN genewise
    ON genewise.locus = all_reciprocal_matches.query_locus
    JOIN taxonomy
    ON taxonomy.species_taxid = assembly_metadata.taxid
);

ALTER TABLE merged ADD COLUMN peptide_length INT;
ALTER TABLE merged ADD COLUMN EVE_length INT;
ALTER TABLE merged ADD COLUMN peptide_coverage INT;
UPDATE merged
	SET peptide_length = CHAR_LENGTH(reconstructed_peptide),
	EVE_length = CAST((SPLIT_PART(SPLIT_PART(EVE_locus, ':', 2), '-', 2)) as INT) - CAST((SPLIT_PART(SPLIT_PART(EVE_locus, ':', 2), '-', 1)) as INT),
	peptide_coverage = CAST((peptide_length * 3) as REAL) / EVE_length * 100;

CREATE TABLE merged_best_hit_rvdb AS (
	SELECT DISTINCT * FROM
		(SELECT *,
		  rank() OVER (PARTITION BY EVE_locus ORDER BY blast_bitscore DESC) ranked
		  FROM (
			SELECT * FROM merged WHERE "database" ILIKE 'rvdb'))
		where ranked = 1) r1;
		
--DROP TABLE assembly_metadata;
--DROP TABLE assembly_statistics;
--DROP TABLE genewise;
--DROP TABLE locus_to_assembly_junction_table;
--DROP TABLE reciprocal_nr;
--DROP TABLE reciprocal_rvdb;
--DROP TABLE taxonomy;