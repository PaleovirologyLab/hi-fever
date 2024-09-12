--Wenyang Lyu, Jose Gabriel Nino Barreat, Emma Harding
--Date: 12/Sep/2024

-- Set schema to use

SET search_path TO postgres;

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

--Create potential output table for downstream analyses, null values coded as N/A

CREATE TABLE potential_output AS
SELECT query_locus, "database",
subject_protein,
COALESCE(SUBSTRING(subject_title FROM 'RecName: Full=(.+?);'),REGEXP_REPLACE(subject_title, '.*\:\s', '')) AS subject_title_clean,
percent_identity, "length", e_value, bitscore,
COALESCE(subject_taxids, 'N/A') AS hit_taxids, COALESCE(superkingdom, 'N/A') AS hit_superkingdom,
COALESCE(kingdom, 'N/A') AS hit_kingdom, COALESCE(phylum, 'N/A') AS hit_phylum,
COALESCE("class", 'N/A') AS "class", COALESCE("order", 'N/A') AS "order",
COALESCE("family", 'N/A') AS "family", COALESCE(genus, 'N/A') AS hit_genus,
COALESCE(species, 'N/A') AS hit_species
FROM all_reciprocal_matches
LEFT JOIN 
ON all_reciprocal_matches.query_locus
ORDER BY query_locus, "database", e_value;

--Create the large merged table

CREATE TABLE merged AS (
SELECT
    all_reciprocal_matches.query_locus as EVE_locus,
    all_reciprocal_matches.percent_identity as BLAST_identity,
    all_reciprocal_matches.e_value as BLAST_evalue,
    all_reciprocal_matches.bitscore as BLAST_bitscore,
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
    genewise.bitscore as genewise_bitscore,
    predicted_ORFs.coverage_of_locus_by_orf,
    predicted_ORFs.orf_seq,
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
    JOIN predicted_ORFs
    ON predicted_ORFs.locus = all_reciprocal_matches.query_locus
    JOIN taxonomy
    ON taxonomy.species_taxid = assembly_metadata.taxid
)

ALTER TABLE merged ADD COLUMN orf_length INT;
UPDATE merged
	SET orf_length = CHAR_LENGTH(orf_seq);
