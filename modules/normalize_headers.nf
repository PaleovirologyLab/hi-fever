process NORMALIZE_ASSEMBLY_HEADERS {

    tag "${assembly}"
    conda 'bioconda::seqkit=2.9.0'

    input:
    path assembly

    output:
    path "${assembly.simpleName}_normalized.fna", emit: normalized_fa
    path "${assembly.simpleName}_header_map.tsv", emit: header_map

    """
    infile="${assembly}"
    if [[ "${assembly}" == *.gz ]]; then
        reader="zcat"
    else
        reader="cat"
    fi

    \$reader "\$infile" | awk '
        BEGIN { OFS="\\t" }
        /^>/ {
            orig = substr(\$0, 2)
            clean = orig
            sub(/ .*/, "", clean)
            sub(/:[0-9]+-[0-9]+$/, "", clean)
            print orig, clean >> "'"${assembly.simpleName}_header_map.tsv"'"
            print ">" clean
            next
        }
        { print }
    ' > "${assembly.simpleName}_normalized.fna"
    """
}
