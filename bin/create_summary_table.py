#!/usr/bin/env python
description = """ 
Combines metadata from sql tables of hi-fever output and classify loci based on hit trends
Usage: python3 /Projects/testing-hi-fever/10genomes_200perfam_noretro/sql/
"""
import modin.pandas as pd
import re
import argparse
import sys
import os

# Unwanted patterns
# un_pat = re.compile(r'isoform.*$|,\spartial|^putative\s', re.IGNORECASE)
# acc_sp = re.compile(r"([A-Z]+_?[0-9]+)(\.[0-9]+)?|\[.*?\]")
# # acc_sp = re.compile(r"([A-Z]+_?[0-9]+)(\.[0-9]+)?|\[.*?\]|endogenous retrovirus")
# rec_name = re.compile(r"RecName: Full=(.+?);")


def clean_subject_title(title):
    """
    Remove unwanted information from hit title
    :param pattern: regular expression of common strings to clean
    :return cleaned string
    """
    pattern = r"([A-Z]+_?[0-9]+)(\.[0-9]+)?|\[.*?\]|LOW QUALITY PROTEIN:|endogenous retrovirus|acc\|[^\|]*\|[^\|]*\|[^\|]*\|[^\|]*\|"
    if pd.isna(title):
        return ""
    cleaned = re.sub(pattern, "", title, flags=re.IGNORECASE)
    return re.sub(r'\s+', ' ', cleaned).strip()

def compute_loci_statistics_nr(data, column="query_locus", order_by="best_bitscore"):
    """
    For each query_locus, return summary metrics
    :param data: modin.pandas.dataframe object, results against nr database
    :param order_by: str, order output
    :return result: modin.pandas.dataframe object, summary counts per loci
    """
    data["eukaryot_hit"] = data["hit_superkingdom"].astype(str).str.contains("eukaryot", case=False)
    data["virus_hit"] = data["hit_superkingdom"].astype(str).str.contains("virus", case=False)
    
    # Best title per locus
    best_titles = data.loc[data.groupby(column)["bitscore"].idxmax()][[column, "cleaned_title", "hit_length"]]
    best_titles = best_titles.rename(columns={
        "cleaned_title": "best_title",
        "hit_length": "best_hit_len"
    })

    # Count of unique species per locus
    n_species = data.groupby(column)["hit_scientific_name"].nunique().reset_index()
    n_species = n_species.rename(columns={"hit_scientific_name": "n_species"})

    # Main stats aggregation
    stats = data.groupby(column).agg(
        n_hits=("query_locus", "count"),
        mean_bitscore=("bitscore", "mean"),
        mean_hit_length=("hit_length", "mean"),
        best_bitscore=("bitscore", "max"),
        common_title=("cleaned_title", lambda x: x.mode().iloc[-1] if not x.mode().empty else ""),
        n_eukaryot_hits=("eukaryot_hit", "sum"),
        n_virus_hits=("virus_hit", "sum")
    ).reset_index()

    # Merge additional info
    result = stats.merge(best_titles, on=column, how="left")
    result = result.merge(n_species, on=column, how="left")

    # Sort by bitscore
    result = result.sort_values(by=order_by, ascending=False)
    return result

def compute_loci_statistics_rvdb(data, column="query_locus", order_by="best_bitscore"):
    """
    For each query_locus, return summary metrics
    :param data: modin.pandas.dataframe object, results against rvdb database 
                                                includes family and kingdom of rvdb hits
    :param order_by: str, order output
    :return result: modin.pandas.dataframe object, summary counts per loci
    """
    data["eukaryot_hit"] = data["hit_superkingdom"].astype(str).str.contains("eukaryot", case=False)
    data["virus_hit"] = data["hit_superkingdom"].astype(str).str.contains("virus", case=False)
    data["retro_hit"] = data["family"].astype(str).str.contains("retro|caulimo|pararna", case=False)

    # Best title per locus
    best_titles = data.loc[data.groupby(column)["bitscore"].idxmax()][[column, "cleaned_title", "hit_length", "family"]]
    best_titles = best_titles.rename(columns={
        "cleaned_title": "best_title",
        "hit_length": "best_hit_len", 
        "family": "best_hit_fam"
    })

    # Count of unique species per locus
    n_species = data.groupby(column)["hit_scientific_name"].nunique().reset_index()
    n_species = n_species.rename(columns={"hit_scientific_name": "n_species"})

    stats = data.groupby(column).agg(
        n_hits=("query_locus", "count"),
        mean_bitscore=("bitscore", "mean"),
        mean_hit_length=("hit_length", "mean"),
        best_bitscore=("bitscore", "max"),
        common_title=("cleaned_title", lambda x: x.mode().iloc[-1] if not x.mode().empty else ""),
        n_eukaryot_hits=("eukaryot_hit", "sum"),
        n_virus_hits=("virus_hit", "sum"),
        n_retros=("retro_hit", "sum"),
        rvdb_fams=("family", lambda x: "|".join(x.dropna().unique())),
        rvdb_kings=("kingdom", lambda x: "|".join(x.dropna().unique()))
    ).reset_index()

    result = stats.merge(best_titles, on=column, how="left")
    result = result.merge(n_species, on=column, how="left")

    # Sort by bitscore
    result = result.sort_values(by=order_by, ascending=False)
    return result

def classify_query_loci(df, rvdb_nr_host_prop =  0.3, retro_euk_prop = 0.7, 
                        rvdb_nr_eve_prop = 0.6, nr_eve_prop = 0.5, rvdb_transposon_prop = 0.8):
    """
    Classification of HI-FEVER results based on:
        - comparisons against nr and rvdb results
        - common labels on host genes
        - common labels on retro-like elements
        - families with high uncertainty
    :param df: modin.pandas.dataframe.DataFrame object, 
                        merged summaries of nr and rvdb results for a single loci
    :param rvdb_nr_host_prop: float [0,1], threshold over which we decide if element is host
                                      i.e, if many nr hits and few rvdb, maybe its a host protein
    :param retro_euk_prop: float [0.1], threshold to define if element is retro from rvdb 
                                        i.e, if proportion of hits to retro plus hits to eukaryot above this threshold
    :param rvdb_nr_eve_prop: float [0.1], threshold to define if element is eve
                                        i.e, if more hits to rvdb than to nr, and not from problematic families
    :param nr_eve_prop: float [0.1], threshold to define if element is eve
                                        i.e, if the majority of hits from nr are viral, then its eve
    :param rvdb_transposon_prop: float [0.1], threshold to define if element is transposon
                                        i.e, if the majority of hits from rvdb are from host, and have certain labels
    :return result: same as input with an additional column for type of element
    """
    def classify(row):
        # === Precompute all needed fields ===
        n_hits_nr = row.get("n_hits_nr", 0)
        n_hits_rvdb = row.get("n_hits_rvdb", 0)
        n_retros = row.get("n_retros", 0)
        n_virus_hits_rvdb = row.get("n_virus_hits_rvdb", 0)
        n_virus_hits_nr = row.get("n_virus_hits_nr", 0)
        n_euk_hits_rvdb = row.get("n_eukaryot_hits_rvdb", 0)
        rvdb_kings = str(row.get("rvdb_kings", "")).lower()
        rvdb_fams = str(row.get("rvdb_fams", "")).lower()

        title_all = " ".join([
            str(row.get("best_title_rvdb", "")),
            str(row.get("common_title_nr", "")),
            str(row.get("common_title_rvdb", ""))
        ]).lower()

        best_title = str(row.get("best_title_rvdb", "")).lower()

        # === Conditions ===
        is_likely_retro = (
            n_hits_rvdb > 0 and
            n_retros > 0 and
            (n_retros + n_euk_hits_rvdb) / n_hits_rvdb > retro_euk_prop and
            re.search(r"gag|pol|env|reverse|pro[- ]?pol|protease|pro|rt[- ]?in", title_all)
        )

        is_host_few_rvdb = (
            n_hits_nr > 0 and
            n_hits_rvdb / n_hits_nr < rvdb_nr_host_prop  # Way more nr than rvdb hits
        )

        is_host_no_rvdb_hits = pd.isna(n_hits_rvdb) or n_hits_rvdb == 0

        is_host_by_name = (
            n_virus_hits_nr == 0 and
            re.search(r"ubiquitin|zinc|nynrin|kinase", title_all)
        )

        is_bamford = re.search(r"bamfordvirae|bamford", rvdb_kings)

        is_transposon = (
            n_hits_rvdb > 0 and
            n_euk_hits_rvdb / n_hits_rvdb > rvdb_transposon_prop and
            re.search(r"\bpol\b(?!yprotein)|transpos|gag|group specific antigen", best_title)
        )

        is_viral_rvdb = (
            n_hits_rvdb > 0 and
            n_virus_hits_rvdb / n_hits_rvdb > rvdb_nr_eve_prop
        )

        is_problem_fam = re.search(r"baculo|poly|fabace|allohe", rvdb_fams)

        is_viral_nr = (
            n_hits_nr > 0 and
            n_virus_hits_nr / n_hits_nr > nr_eve_prop
        )

        # === Decision tree ===
        if is_likely_retro:
            return "likely-retro"
        elif is_host_few_rvdb:
            return "likely-host-protein"
        elif is_host_no_rvdb_hits:
            return "likely-host-protein"
        elif is_host_by_name:
            return "likely-host-protein"
        elif is_bamford:
            return "likely-bamfordviridae"
        elif is_transposon:
            return "likely-transposon"
        elif is_viral_rvdb:
            return "likely-eve" if not is_problem_fam else "uncertain"
        elif is_viral_nr:
            return "likely-eve"
        else:
            return "uncertain"

    # Apply function row-wise (use df._to_pandas() if in Modin and failing)
    result = df.copy()
    result["element_type"] = result.apply(classify, axis=1)
    return result


def filter_rvdb_by_taxonomy(rvdb_table, uniq_tax, pattern):
    temp = uniq_tax['family'].str.contains(pattern, flags=re.IGNORECASE, regex=True, na=False)
    filtered_tax = uniq_tax[temp]
    rvdb_res_filtered = pd.merge(rvdb_table, filtered_tax, left_on="taxid", right_on="taxid", suffixes=("", "_rvdb"), how="inner")
    return (rvdb_res_filtered)

def resolve_path_or_error(base, filename):
    path = os.path.join(base, filename)
    if not os.path.exists(path):
        raise FileNotFoundError(f"Expected file '{filename}' not found in folder '{base}'")
    return path

def parse_arguments():
    parser = argparse.ArgumentParser(description=description)

    # Optional folder input
    parser.add_argument("--folder", type=str, help="Path to folder with all input files")

    # Direct file inputs
    parser.add_argument("--reciprocal_nr", type=str)
    parser.add_argument("--reciprocal_rvdb", type=str)
    parser.add_argument("--taxonomy", type=str)
    parser.add_argument("--assembly_map", type=str)
    parser.add_argument("--assembly_metadata", type=str)
    parser.add_argument("--genewise", type=str)

    # Options
    parser.add_argument("-al", "--aa_len", type=int, default=1, help="Min peptide length")
    parser.add_argument("-l", "--label", type=str, default='summary', help="Label for output")
    parser.add_argument("-sf", "--sequence_file", action='store_false', help="Disable FASTA output")
    parser.add_argument("-fb", "--filter_by", type=str, default=False)
    parser.add_argument("-ov", "--old_version", action='store_true', help="All version had different notation for hits taxonomy")

    args = parser.parse_args()

    # Infer file paths from folder if necessary
    if args.folder:
        args.reciprocal_nr = args.reciprocal_nr or resolve_path_or_error(args.folder, "reciprocal-nr-matches.dmnd.tsv")
        args.reciprocal_rvdb = args.reciprocal_rvdb or resolve_path_or_error(args.folder, "reciprocal-rvdb-matches.dmnd.tsv")
        args.taxonomy = args.taxonomy or resolve_path_or_error(args.folder, "hits_taxonomy.tsv")
        args.assembly_map = args.assembly_map or resolve_path_or_error(args.folder, "locus_assembly_map.tsv")
        args.assembly_metadata = args.assembly_metadata or resolve_path_or_error(args.folder, "assembly_metadata.tsv")
        args.genewise = args.genewise or resolve_path_or_error(args.folder, "genewise.tsv")
    else:
        args.folder = os.path.dirname(args.reciprocal_nr)

    # Final sanity check
    missing = [k for k in ['reciprocal_nr', 'reciprocal_rvdb', 'taxonomy', 'assembly_map', 'assembly_metadata', 'genewise']
               if getattr(args, k) is None]
    if missing:
        parser.error(f"Missing required files: {', '.join(missing)}. Provide --folder or set them individually.")

    return args

def main():
    args = parse_arguments()

    # === Load reciprocal data
    reciprocal_columns = [
        'query_locus', 'subject_protein', 'percent_identity', 'hit_length', 'mismatches', 'gapopens',
        'query_start', 'query_end', 'subject_start', 'subject_end', 'e_value', 'bitscore',
        'hit_taxid', 'hit_scientific_name', 'hit_superkingdom', 'hit_kingdom',
        'hit_phylum', 'hit_title'
    ]

    d_nr = pd.read_csv(args.reciprocal_nr, sep="\t", header=None, names=reciprocal_columns)
    d_rvdb = pd.read_csv(args.reciprocal_rvdb, sep="\t", header=None, names=reciprocal_columns)

    d_nr["cleaned_title"] = d_nr["hit_title"].apply(clean_subject_title)
    d_rvdb["cleaned_title"] = d_rvdb["hit_title"].apply(clean_subject_title)

    tax = pd.read_csv(args.taxonomy, sep="\t", header=None,
                      names=['taxid', 'superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
    
    # Coerce bad entries to NaN, then drop or fill
    d_rvdb["hit_taxid_int"] = pd.to_numeric(d_rvdb["hit_taxid"], errors="coerce").fillna(-1).astype(int)
    tax["taxid_int"] = pd.to_numeric(tax["taxid"], errors="coerce").fillna(-1).astype(int)

    d_rvdb_t = d_rvdb.merge(tax, left_on="hit_taxid_int", right_on="taxid_int", how="left")
    # d_rvdb_t = d_rvdb.merge(tax, left_on="hit_taxid", right_on="taxid", how="left")
    
    # Old version of HI-FEVER had a different convention for storing taxonomical information 
    if args.old_version:
        # Ensure both are strings (safe, avoids casting errors)
        d_nr["hit_taxid_str"] = d_nr["hit_taxid"].astype(str)
        tax["taxid_str"] = tax["taxid"].astype(str)

        # Create a lookup dictionary from tax table
        tax_lookup = tax.set_index("taxid_str")["superkingdom"].to_dict()

        # Map superkingdom to d_nr based on hit_taxid_str
        d_nr["hit_superkingdom"] = d_nr["hit_taxid_str"].map(tax_lookup)

    loci_stats_nr = compute_loci_statistics_nr(d_nr)
    d_rvdb_t["family"] = d_rvdb_t["family"].astype(str)
    d_rvdb_t = d_rvdb_t.sort_values(by="family", kind="mergesort", ignore_index=True)
    loci_stats_rvdb = compute_loci_statistics_rvdb(d_rvdb_t)

    both_dbs = loci_stats_nr.merge(loci_stats_rvdb, on="query_locus", how="outer", suffixes=("_nr", "_rvdb"))
    res = classify_query_loci(both_dbs)

    # === Assembly metadata
    a_map = pd.read_csv(args.assembly_map, sep='\t', header=None, names=['query_locus', 'assembly'])
    a_map["assembly"] = a_map["assembly"].astype(str).str.extract(r"\b([A-Z]{3}_[0-9]+(?:\.[0-9]+)?)")

    metadata = pd.read_csv(args.assembly_metadata, sep='\t', usecols=[0, 1], header=None, names=['hostName', 'assembly'])
    host = pd.merge(a_map, metadata, on="assembly", how="inner")

    # === Genewise data
    genewise_res = pd.read_csv(args.genewise, sep='\t', usecols=[4, 10, 11, 12, 13, 14], header=None,
                               names=['query_locus', 'cdna_seq', 'peptide_seq', 'introns', 'frameshifts', 'stop_codons'])
    genewise_res['query_locus'] = genewise_res['query_locus'].astype('string')
    genewise_res['peptide_len'] = genewise_res['peptide_seq'].str.len()
    filtered_gw = genewise_res[genewise_res.peptide_len > args.aa_len]

    # === Merge all
    merged = pd.merge(filtered_gw, host, on="query_locus", suffixes=("_gw", "_host"), how="left")
    final = pd.merge(res, merged, on="query_locus", how="left")
    final["peptideID"] = ">" + final[['query_locus', 'hostName', 'best_title_rvdb']].fillna("").astype(str).agg('|'.join, axis=1)

    output_prefix = os.path.join(args.folder, args.label)

    merged_no_seq = final.drop(['peptide_seq', 'peptideID', 'cdna_seq'], axis=1)
    merged_no_seq.to_csv(f'{output_prefix}_classified.tsv', index=False, sep='\t')

    # === FASTA output
    if args.sequence_file:
        # Get only reconstructed sequences
        final = final[final['peptide_len'].notna() & (final['peptide_len'] > 0)]

        aa_fasta = final.loc[:, ["peptideID", "peptide_seq"]]
        aa_fasta.to_csv(f'{output_prefix}_aa_genewise.fasta', index=False, sep='\n', header=None)

        nt_fasta = final.loc[:, ["peptideID", "cdna_seq"]]
        nt_fasta.to_csv(f'{output_prefix}_cdna_genewise.fasta', index=False, sep='\n', header=None)

if __name__ == "__main__":
    main()