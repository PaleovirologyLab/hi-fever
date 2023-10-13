#!/bin/bash -ue
export WISECONFIGDIR="$CONDA_PREFIX/share/wise2/wisecfg"
mkdir wise_tmp wise_tmp_nt wise_tmp_nt2 wise_tmp_aa

# Generate table of original best proteins per locus

cut -f5,8 pre_reciprocal_genewise.txt |     sort -k1,1 >     wise_tmp/original_pairs

# Generate table of reciprocal best proteins per locus

sort -k1,1 -k12,12nr reciprocal-matches.dmnd.tsv |     sort -u -k1,1 |     cut -f1-2 >     wise_tmp/reciprocal_pairs

# Find loci whose best hit has changed & get nuc + protein identifiers + genomic coords & file to intersect with context coords

comm -13 wise_tmp/original_pairs wise_tmp/reciprocal_pairs |     tee >(cut -f2 | sort | uniq > wise_tmp/protein_accessions) |     tee >(cut -f1 | sed 's/:/	/; s/-/	/' > wise_tmp/genomic_coords) |     tee >(sed 's/:/	/; s/-/	/' | sort -k1,1 -k2,2n > wise_tmp/intersection_bed) >     wise_tmp/matched_pairs

# Initialise nextflow file path (due to nextflow bug?)

touch best_reciprocals.fasta

# Extract new best proteins into individual files

seqtk subseq best_reciprocals.fasta wise_tmp/protein_accessions |     sed 's/ .*//' >     wise_tmp/temp.fa

awk '/^>/ { if (name) close(name); name="wise_tmp_aa/" substr($0,2); print > name; next } { print >> name }' wise_tmp/temp.fa

rm "$(readlink -f best_reciprocals.fasta)"

# Get all strict FASTAs together and extract relevant ones

cat GCA_023065795.1_ASM2306579v1_genomic.fna_strict.fasta GCA_000601435.1_P_lyc_v1_genomic.fna_strict.fasta GCA_026225775.1_INECOL_Bpal_1.0_genomic.fna_strict.fasta GCA_013400295.1_ASM1340029v1_genomic.fna_strict.fasta GCA_025773905.1_ASM2577390v1_genomic.fna_strict.fasta GCA_004802705.1_ASM480270v1_genomic.fna_strict.fasta GCA_016432905.1_ASM1643290v1_genomic.fna_strict.fasta GCA_000006255.1_Postia_placenta_V1.0_genomic.fna_strict.fasta GCA_014048405.2_ASM1404840v2_genomic.fna_strict.fasta GCA_005959805.1_ASM595980v1_genomic.fna_strict.fasta GCA_004348175.1_PNOVO_Noble_KM_genomic.fna_strict.fasta GCA_000467735.1_ASM46773v1_genomic.fna_strict.fasta GCA_027569095.1_ASM2756909v1_genomic.fna_strict.fasta GCA_003243085.1_ASM324308v1_genomic.fna_strict.fasta GCA_000226015.1_MneLei_Aug2011_genomic.fna_strict.fasta GCA_014885075.1_ASM1488507v1_genomic.fna_strict.fasta GCA_022814215.1_ASM2281421v1_genomic.fna_strict.fasta GCA_000330505.1_EIA2_v2_genomic.fna_strict.fasta GCA_900406335.1_ASM90040633v1_genomic.fna_strict.fasta GCA_002846955.1_A_aquasalis_v1.0_genomic.fna_strict.fasta GCA_002076135.1_ASM207613v1_genomic.fna_strict.fasta GCA_027577635.1_UCR_CoemRSA475_1.0_genomic.fna_strict.fasta GCA_018691285.1_UNAM_Csor_1.0_genomic.fna_strict.fasta GCA_000181375.1_ASM18137v1_genomic.fna_strict.fasta GCA_002891735.1_TetSoc1_genomic.fna_strict.fasta GCA_030078075.1_ASM3007807v1_genomic.fna_strict.fasta GCA_000164785.2_C_hoffmanni-2.0.1_genomic.fna_strict.fasta GCA_015342145.1_ASM1534214v1_genomic.fna_strict.fasta GCA_009602735.1_ASM960273v1_genomic.fna_strict.fasta GCA_024256425.2_CGAR_prim_01v2_genomic.fna_strict.fasta GCA_002938485.2_Soryzae_2.0_genomic.fna_strict.fasta GCA_002018215.1_CTBE_SP803280_v1.0_genomic.fna_strict.fasta GCA_000317765.2_CHIR_2.0_genomic.fna_strict.fasta GCA_022539395.1_myiCay2020_genomic.fna_strict.fasta GCA_011064425.1_Rrattus_CSIRO_v1_genomic.fna_strict.fasta GCA_003571905.1_ASM357190v1_genomic.fna_strict.fasta GCA_027569235.1_ASM2756923v1_genomic.fna_strict.fasta GCA_026283905.1_ASM2628390v1_genomic.fna_strict.fasta GCA_000152225.2_Pcap_2.0_genomic.fna_strict.fasta GCA_010312235.1_ButJap1.0_genomic.fna_strict.fasta GCA_027575695.1_UCR_CoemRSA638_1.0_genomic.fna_strict.fasta GCA_018249875.1_SRR7174573_genomic.fna_strict.fasta GCF_000002285.5_Dog10K_Boxer_Tasha_genomic.fna_strict.fasta GCA_002966915.1_ASM296691v1_genomic.fna_strict.fasta GCA_023638675.1_C.penn_1.0_hap1_genomic.fna_strict.fasta GCA_013401355.1_ASM1340135v1_genomic.fna_strict.fasta GCA_001599015.1_JCM_9580_assembly_v001_genomic.fna_strict.fasta GCA_001239805.1_ASM123980v1_genomic.fna_strict.fasta GCF_000001635.27_GRCm39_genomic.fna_strict.fasta GCF_000296755.1_EriEur2.0_genomic.fna_strict.fasta GCA_002116315.1_ASM211631v1_genomic.fna_strict.fasta GCA_022495145.1_Xylbam139988_1_genomic.fna_strict.fasta GCA_022651745.1_ASM2265174v1_genomic.fna_strict.fasta GCA_025766325.1_Benpoi1_genomic.fna_strict.fasta GCF_000313985.2_ASM31398v2_genomic.fna_strict.fasta GCA_022817305.1_ASM2281730v1_genomic.fna_strict.fasta GCA_013396695.1_ASM1339669v1_genomic.fna_strict.fasta GCA_003013575.1_ASM301357v1_genomic.fna_strict.fasta GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.fna_strict.fasta GCA_903813345.2_Owenia_chromosome_genomic.fna_strict.fasta GCA_911387705.1_ihAelAcum1.1_alternate_haplotype_genomic.fna_strict.fasta GCF_000005575.2_AgamP3_genomic.fna_strict.fasta GCA_018119265.1_PriViv1.0_genomic.fna_strict.fasta GCA_017591425.1_NIBS_Ocur_1.0_genomic.fna_strict.fasta GCF_011100555.1_mCalJa1.2.pat.X_genomic.fna_strict.fasta GCA_023781915.1_ASM2378191v1_genomic.fna_strict.fasta GCA_948208565.1_Hc-NTBG20172064-DRAFT-NextGenCassava-1.0_genomic.fna_strict.fasta GCF_000147115.1_Myoluc2.0_genomic.fna_strict.fasta GCF_000003025.6_Sscrofa11.1_genomic.fna_strict.fasta GCF_000001405.40_GRCh38.p14_genomic.fna_strict.fasta GCF_000151845.1_Pvam_2.0_genomic.fna_strict.fasta GCF_000001905.1_Loxafr3.0_genomic.fna_strict.fasta GCF_000151885.1_Dord_2.0_genomic.fna_strict.fasta GCF_000151735.1_Cavpor3.0_genomic.fna_strict.fasta GCF_000164845.4_VicPac3.2_genomic.fna_strict.fasta GCF_000165445.2_Mmur_3.0_genomic.fna_strict.fasta GCF_016920785.2_ASM1692078v2_genomic.fna_strict.fasta GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna_strict.fasta GCF_000181295.1_OtoGar3_genomic.fna_strict.fasta GCF_000002295.2_MonDom5_genomic.fna_strict.fasta GCF_002863925.1_EquCab3.0_genomic.fna_strict.fasta GCF_002263795.2_ARS-UCD1.3_genomic.fna_strict.fasta GCF_000164805.1_Tarsius_syrichta-2.0.1_genomic.fna_strict.fasta GCF_004115215.2_mOrnAna1.pri.v4_genomic.fna_strict.fasta GCF_002007445.2_ASM200744v3_genomic.fna_strict.fasta GCF_015227675.2_mRatBN7.2_genomic.fna_strict.fasta GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.fna_strict.fasta GCF_000208655.3_Dasnov3.2_genomic.fna_strict.fasta GCA_028372415.1_mMacEug1.pri_genomic.fna_strict.fasta GCF_003957565.2_bTaeGut1.4.pri_genomic.fna_strict.fasta GCF_014633375.1_OchPri4.0_genomic.fna_strict.fasta GCF_003339765.1_Mmul_10_genomic.fna_strict.fasta GCF_002204515.2_AaegL5.0_genomic.fna_strict.fasta GCF_011762595.1_mTurTru1.mat.Y_genomic.fna_strict.fasta GCF_016881025.1_HiC_Itri_2_genomic.fna_strict.fasta GCF_009806435.1_UM_NZW_1.0_genomic.fna_strict.fasta GCF_027595985.1_mSorAra2.pri_genomic.fna_strict.fasta GCF_028858775.1_NHGRI_mPanTro3-v1.1-hic.freeze_pri_genomic.fna_strict.fasta GCF_029281585.1_NHGRI_mGorGor1-v1.1-0.2.freeze_pri_genomic.fna_strict.fasta GCF_028885625.1_NHGRI_mPonPyg2-v1.1-hic.freeze_pri_genomic.fna_strict.fasta |     sed 's/ .*//' >     wise_tmp/temp.fa

awk '/^>/ { if (name) close(name); name="wise_tmp_nt/" substr($0,2); print > name; next } { print >> name }' wise_tmp/temp.fa

# GeneWise operations, strict FASTA

while read line
    do
        query=$(echo $line | cut -f2 -d ' ')
        target=$(echo $line | cut -f1 -d ' ')
        genewise wise_tmp_aa/$query wise_tmp_nt/$target -both -matrix "BLOSUM62".bla -sum -pep -cdna -divide DIVIDE_STRING -silent |             grep -v ">\|Bits   Query" |             awk '/^[-0-9]/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' |             sed 's/DIVIDE_STRING/	/g' |             tr -s '	' |             tr -s ' ' '	' |             sort -k5,5 -k1,1nr |             sort -u -k5,5 >>             wise_tmp/genewise_strict
    done < wise_tmp/matched_pairs

# Back-calculate genomic coords from genewise

paste wise_tmp/genewise_strict wise_tmp/genomic_coords |     tr -s '	' |     awk 'BEGIN{OFS="	"} {if ($6 < $7) print $12, $13+$6-1, $13+$7-1, "+", $5, "strict_PR", $1, $2, $3, $4, $8, $9, $10, $11; else print $12, $13+$7-1, $13+$6-1, "-", $5, "strict_PR", $1, $2, $3, $4, $8, $9, $10, $11}' |     sort -k1,1 -k2,2n >     wise_tmp/output1

# Get all context FASTAs together

cat GCA_023065795.1_ASM2306579v1_genomic.fna_context.fasta GCA_000601435.1_P_lyc_v1_genomic.fna_context.fasta GCA_026225775.1_INECOL_Bpal_1.0_genomic.fna_context.fasta GCA_013400295.1_ASM1340029v1_genomic.fna_context.fasta GCA_025773905.1_ASM2577390v1_genomic.fna_context.fasta GCA_004802705.1_ASM480270v1_genomic.fna_context.fasta GCA_016432905.1_ASM1643290v1_genomic.fna_context.fasta GCA_000006255.1_Postia_placenta_V1.0_genomic.fna_context.fasta GCA_014048405.2_ASM1404840v2_genomic.fna_context.fasta GCA_005959805.1_ASM595980v1_genomic.fna_context.fasta GCA_004348175.1_PNOVO_Noble_KM_genomic.fna_context.fasta GCA_000467735.1_ASM46773v1_genomic.fna_context.fasta GCA_027569095.1_ASM2756909v1_genomic.fna_context.fasta GCA_003243085.1_ASM324308v1_genomic.fna_context.fasta GCA_000226015.1_MneLei_Aug2011_genomic.fna_context.fasta GCA_014885075.1_ASM1488507v1_genomic.fna_context.fasta GCA_022814215.1_ASM2281421v1_genomic.fna_context.fasta GCA_000330505.1_EIA2_v2_genomic.fna_context.fasta GCA_900406335.1_ASM90040633v1_genomic.fna_context.fasta GCA_002846955.1_A_aquasalis_v1.0_genomic.fna_context.fasta GCA_002076135.1_ASM207613v1_genomic.fna_context.fasta GCA_027577635.1_UCR_CoemRSA475_1.0_genomic.fna_context.fasta GCA_018691285.1_UNAM_Csor_1.0_genomic.fna_context.fasta GCA_000181375.1_ASM18137v1_genomic.fna_context.fasta GCA_002891735.1_TetSoc1_genomic.fna_context.fasta GCA_030078075.1_ASM3007807v1_genomic.fna_context.fasta GCA_000164785.2_C_hoffmanni-2.0.1_genomic.fna_context.fasta GCA_015342145.1_ASM1534214v1_genomic.fna_context.fasta GCA_009602735.1_ASM960273v1_genomic.fna_context.fasta GCA_024256425.2_CGAR_prim_01v2_genomic.fna_context.fasta GCA_002938485.2_Soryzae_2.0_genomic.fna_context.fasta GCA_002018215.1_CTBE_SP803280_v1.0_genomic.fna_context.fasta GCA_000317765.2_CHIR_2.0_genomic.fna_context.fasta GCA_022539395.1_myiCay2020_genomic.fna_context.fasta GCA_011064425.1_Rrattus_CSIRO_v1_genomic.fna_context.fasta GCA_003571905.1_ASM357190v1_genomic.fna_context.fasta GCA_027569235.1_ASM2756923v1_genomic.fna_context.fasta GCA_026283905.1_ASM2628390v1_genomic.fna_context.fasta GCA_000152225.2_Pcap_2.0_genomic.fna_context.fasta GCA_010312235.1_ButJap1.0_genomic.fna_context.fasta GCA_027575695.1_UCR_CoemRSA638_1.0_genomic.fna_context.fasta GCA_018249875.1_SRR7174573_genomic.fna_context.fasta GCF_000002285.5_Dog10K_Boxer_Tasha_genomic.fna_context.fasta GCA_002966915.1_ASM296691v1_genomic.fna_context.fasta GCA_023638675.1_C.penn_1.0_hap1_genomic.fna_context.fasta GCA_013401355.1_ASM1340135v1_genomic.fna_context.fasta GCA_001599015.1_JCM_9580_assembly_v001_genomic.fna_context.fasta GCA_001239805.1_ASM123980v1_genomic.fna_context.fasta GCF_000001635.27_GRCm39_genomic.fna_context.fasta GCF_000296755.1_EriEur2.0_genomic.fna_context.fasta GCA_002116315.1_ASM211631v1_genomic.fna_context.fasta GCA_022495145.1_Xylbam139988_1_genomic.fna_context.fasta GCA_022651745.1_ASM2265174v1_genomic.fna_context.fasta GCA_025766325.1_Benpoi1_genomic.fna_context.fasta GCF_000313985.2_ASM31398v2_genomic.fna_context.fasta GCA_022817305.1_ASM2281730v1_genomic.fna_context.fasta GCA_013396695.1_ASM1339669v1_genomic.fna_context.fasta GCA_003013575.1_ASM301357v1_genomic.fna_context.fasta GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.fna_context.fasta GCA_903813345.2_Owenia_chromosome_genomic.fna_context.fasta GCA_911387705.1_ihAelAcum1.1_alternate_haplotype_genomic.fna_context.fasta GCF_000005575.2_AgamP3_genomic.fna_context.fasta GCA_018119265.1_PriViv1.0_genomic.fna_context.fasta GCA_017591425.1_NIBS_Ocur_1.0_genomic.fna_context.fasta GCF_011100555.1_mCalJa1.2.pat.X_genomic.fna_context.fasta GCA_023781915.1_ASM2378191v1_genomic.fna_context.fasta GCA_948208565.1_Hc-NTBG20172064-DRAFT-NextGenCassava-1.0_genomic.fna_context.fasta GCF_000147115.1_Myoluc2.0_genomic.fna_context.fasta GCF_000003025.6_Sscrofa11.1_genomic.fna_context.fasta GCF_000001405.40_GRCh38.p14_genomic.fna_context.fasta GCF_000151845.1_Pvam_2.0_genomic.fna_context.fasta GCF_000001905.1_Loxafr3.0_genomic.fna_context.fasta GCF_000151885.1_Dord_2.0_genomic.fna_context.fasta GCF_000151735.1_Cavpor3.0_genomic.fna_context.fasta GCF_000164845.4_VicPac3.2_genomic.fna_context.fasta GCF_000165445.2_Mmur_3.0_genomic.fna_context.fasta GCF_016920785.2_ASM1692078v2_genomic.fna_context.fasta GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna_context.fasta GCF_000181295.1_OtoGar3_genomic.fna_context.fasta GCF_000002295.2_MonDom5_genomic.fna_context.fasta GCF_002863925.1_EquCab3.0_genomic.fna_context.fasta GCF_002263795.2_ARS-UCD1.3_genomic.fna_context.fasta GCF_000164805.1_Tarsius_syrichta-2.0.1_genomic.fna_context.fasta GCF_004115215.2_mOrnAna1.pri.v4_genomic.fna_context.fasta GCF_002007445.2_ASM200744v3_genomic.fna_context.fasta GCF_015227675.2_mRatBN7.2_genomic.fna_context.fasta GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.fna_context.fasta GCF_000208655.3_Dasnov3.2_genomic.fna_context.fasta GCA_028372415.1_mMacEug1.pri_genomic.fna_context.fasta GCF_003957565.2_bTaeGut1.4.pri_genomic.fna_context.fasta GCF_014633375.1_OchPri4.0_genomic.fna_context.fasta GCF_003339765.1_Mmul_10_genomic.fna_context.fasta GCF_002204515.2_AaegL5.0_genomic.fna_context.fasta GCF_011762595.1_mTurTru1.mat.Y_genomic.fna_context.fasta GCF_016881025.1_HiC_Itri_2_genomic.fna_context.fasta GCF_009806435.1_UM_NZW_1.0_genomic.fna_context.fasta GCF_027595985.1_mSorAra2.pri_genomic.fna_context.fasta GCF_028858775.1_NHGRI_mPanTro3-v1.1-hic.freeze_pri_genomic.fna_context.fasta GCF_029281585.1_NHGRI_mGorGor1-v1.1-0.2.freeze_pri_genomic.fna_context.fasta GCF_028885625.1_NHGRI_mPonPyg2-v1.1-hic.freeze_pri_genomic.fna_context.fasta |     sed 's/ .*//' >     wise_tmp/temp.fa

#  Convert context ranges to BED for intersection with those strict regions where the best protein changed

grep ">" wise_tmp/temp.fa |     sed 's/>//; s/:/	/; s/-/	/' |     sort -k1,1 -k2,2n >     wise_tmp/all_context_bed

# Intersect ranges, keep various outputs from the context ranges we need

bedtools intersect -a wise_tmp/intersection_bed -b wise_tmp/all_context_bed -wb -f 1 -sorted |     tee >(cut -f5-7 > wise_tmp/genomic_coords) |     awk 'BEGIN{OFS="	"} {print $5":"$6"-"$7, $4}' >     wise_tmp/matched_pairs

# Extract those sequences

awk '/^>/ { if (name) close(name); name="wise_tmp_nt2/" substr($0,2); print > name; next } { print >> name }' wise_tmp/temp.fa

# GeneWise operations, context FASTA

while read line
    do
        query=$(echo $line | cut -f2 -d ' ')
        target=$(echo $line | cut -f1 -d ' ')
        genewise wise_tmp_aa/$query wise_tmp_nt2/$target -both -matrix "BLOSUM62".bla -sum -pep -cdna -divide DIVIDE_STRING -silent |             grep -v ">\|Bits   Query" |             awk '/^[-0-9]/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' |             sed 's/DIVIDE_STRING/	/g' |             tr -s '	' |             tr -s ' ' '	' |             sort -k5,5 -k1,1nr |             sort -u -k5,5 >>             wise_tmp/genewise_context
    done < wise_tmp/matched_pairs

# Back-calculate genomic coords from genewise & remove redundancy (sites found in two context FASTAs)

paste wise_tmp/genewise_context wise_tmp/genomic_coords |     tr -s '	' |     awk 'BEGIN{OFS="	"} {if ($6 < $7) print $12, $13+$6-1, $13+$7-1, "+", $5, "context_PR", $1, $2, $3, $4, $8, $9, $10, $11; else print $12, $13+$7-1, $13+$6-1, "-", $5, "context_PR", $1, $2, $3, $4, $8, $9, $10, $11}' |     sort -u -k1,1 -k2,2n -k3,3nr >     wise_tmp/output2

# Intersect and concatenate results (keep strict if not covered by context, otherwise keep context)

# First report strict predictions not encompassed by a context prediction

bedtools intersect -v -a wise_tmp/output1 -b wise_tmp/output2 -f 1 -wa > wise_tmp/merged_results

# Second find strict predictions encompassed by context predictions, report latter

bedtools intersect -a wise_tmp/output1 -b wise_tmp/output2 -f 1 -wb |     awk 'BEGIN{OFS="	"} {print $15, $16, $17, $18, $5, $20, $21, $22, $23, $24, $25, $26, $27, $28}' >>     wise_tmp/merged_results

# Post-processing of in-frame STOPs

stop_convert_and_count.py --task remove --file wise_tmp/merged_results > post_reciprocal_genewise.txt

# Merge both genewise outputs - for each locus keep the longest individual prediction

cat pre_reciprocal_genewise.txt post_reciprocal_genewise.txt |     sort -k1,1 -k2,2n -k3,3nr |     sort -u -k5,5 >     genewise.txt

# Cleanup

rm -r wise_tmp*
