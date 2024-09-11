--Wenyang Lyu, Jose Gabriel Nino Barreat
--Date: 12/Dec/2023

-- Set schema to use

SET search_path TO hifever_schema;

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

--Subset loci names and model labels from the domain hits

CREATE TABLE domains AS
SELECT locus,
ARRAY_TO_STRING(ARRAY_AGG(domains_i_evalue_best_single), ', ', 'N/A') AS domains_i_evalue_best,
ARRAY_TO_STRING(ARRAY_AGG(domains_bitscore_best_single), ', ', 'N/A') AS domains_bitscore_best,
ARRAY_TO_STRING(ARRAY_AGG(model_name), ', ', 'N/A') AS model_names,
ARRAY_TO_STRING(ARRAY_AGG(model_description), ', ', 'N/A') AS model_descriptions,
ARRAY_TO_STRING(ARRAY_AGG(model_acc), ', ', 'N/A') AS model_accs
FROM best_forward_hits_phmm
GROUP BY locus;

--Create potential output table for downstream analyses, null values coded as N/A

CREATE TABLE potential_output AS
SELECT query_locus, "database",
subject_protein,
--subject_title,
REGEXP_REPLACE(COALESCE(SUBSTRING(subject_title FROM 'RecName: Full=(.+?);'),REGEXP_REPLACE(subject_title, '.*\:\s', '')), 'isoform.*$|,\spartial|^putative\s', '') AS subject_title_clean,
domains_i_evalue_best, 
domains_bitscore_best, model_names, model_descriptions, model_accs, percent_identity, "length", mismatches,
gapopens, e_value, bitscore,
COALESCE(subject_taxids, 'N/A') AS subject_taxids, COALESCE(superkingdom, 'N/A') AS superkingdom,
COALESCE(kingdom, 'N/A') AS kingdom, COALESCE(phylum, 'N/A') AS phylum,
COALESCE("class", 'N/A') AS "class", COALESCE("order", 'N/A') AS "order",
COALESCE("family", 'N/A') AS "family", COALESCE(genus, 'N/A') AS genus,
COALESCE(species, 'N/A') AS species
FROM all_reciprocal_matches
LEFT JOIN domains
ON all_reciprocal_matches.query_locus = domains.locus
ORDER BY query_locus, "database", e_value;

SELECT * FROM potential_output
WHERE "family" LIKE 'Filoviridae';

DROP TABLE all_reciprocal_matches;
DROP TABLE potential_output;
DROP TABLE domains;

SELECT * FROM all_reciprocal_matches;

SELECT * FROM domains;

SELECT * FROM best_forward_hits_phmm
WHERE locus LIKE 'AAPY01115927%';

SELECT * FROM best_forward_hits_phmm;
SELECT * FROM reciprocal_nr;
