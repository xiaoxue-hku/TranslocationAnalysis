#!/bin/sh

# reference genome hg19 

# genomic region files downloaded from UCSC Table Browser 

# gene expression

declare -a arr=("SARC" "BRCA" "COAD" "OV" "LUSC" "STAD" "UCEC" "HNSC" "LIHC" "SKCM")
for cancer in "${arr[@]}"
do
	wget "https://gdc.xenahubs.net/download/TCGA-$cancer.htseq_fpkm-uq.tsv.gz" && gunzip TCGA-$cancer.htseq_fpkm-uq.tsv.gz
done


# HEPG2 RepliSeq data 

intersectBed -a hg19.knownCanonical.autosome.bed.txt -b HEPG2_RepliSeq.bed -wo | sortBed -i - | uniq -  > hg19_HEPG2.rep.bed.txt && datamash -s -g 4 mean 10 < hg19_HEPG2.rep.bed.txt > hg19_HEPG2.rep.tsv

awk -F '\t' '{print $0 "\t" "In"NR "\t" "." "\t" "."}' hg19.intergenic-2kb.autosome.bed.txt | sortBed -i - | uniq - > hg19_numbered.intergenic-2kb.autosome.bed.txt
intersectBed -a hg19_numbered.intergenic-2kb.autosome.bed.txt -b HEPG2_RepliSeq.bed -wo | sortBed -i - | uniq -  > hg19_HEPG2_intergenic.rep.bed.txt && datamash -s -g 4 mean 10 < hg19_HEPG2_intergenic.rep.bed.txt > hg19_HEPG2_intergenic.rep.tsv

bedtools makewindows -g hg19.chrom.sizes -w 1000000 | sed '/chrX\|chrY\|chrUn\|chrM\|random\|hap/d' - | sortBed -i - | uniq - | awk -F '\t' '{print $0 "\t" "oneMb"NR "\t" "." "\t" "."}' - > hg19_numbered.1Mb.autosome.bed.txt
intersectBed -a hg19_numbered.1Mb.autosome.bed.txt -b HEPG2_RepliSeq.bed -wo | sortBed -i - | uniq -  > hg19_HEPG2_1Mb.rep.bed.txt && datamash -s -g 4 mean 10 < hg19_HEPG2_1Mb.rep.bed.txt > hg19_HEPG2_1Mb.rep.tsv
cut -f3,4 hg19_numbered.1Mb.autosome.bed.txt | sort - | uniq  - > convertID_1Mb.tsv


# structural somatic mutation from raw PCAWG (ICGC+TCGA) database
wget https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_sv/final_consensus_sv_bedpe_passonly.icgc.public.tgz && wget https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_sv/final_consensus_sv_bedpe_passonly.tcga.public.tgz

cat icgc_allTRA.bed tcga_allTRA.bed | sed '/chrX\|chrY/d' - | sortBed -i - > rawPCAWG_allTRA.bed	# combine extracted translocations from both icgc and tcga for following analysis

awk -F '\t' 'NR==FNR{xena[$1]=$4;next} $2 in xena {print $0 "\t" xena[$2]}' PCAWG_sub_signatures_in_samples_beta2.20170320.donor pcawg_sample_info.txt > pcawg_sample_Sig3.txt

awk -F '\t' 'NR==FNR{raw[$1]=$2 "\t" $3;next} $4 in raw {print $1 "\t" $2 "\t" $3 "\t" "." "\t" "." "\t" "TRA" "\t" raw[$4]}' pcawg_sample_info.txt rawPCAWG_allTRA.bed | cut -f1 -d '-' - | sed 's/BOCA/SARC/g' - | sortBed -i - | uniq - > rawPCAWG_allTRA_ann.bed

declare -a arr=("SARC" "BRCA" "GBM" "OV" "LUSC" "STAD" "UCEC" "HNSC" "LIHC" "SKCM")
for cancer in "${arr[@]}"
do
	grep -F "$cancer" rawPCAWG_allTRA_ann.bed | sortBed -i - | uniq - > "$cancer"_TRA.bed
done


# mutational signature

wget https://pcawg.xenahubs.net/download/PCAWG_sub_signatures_in_samples_beta2.20170320.donor


# PCAWG project code

wget https://pcawg.xenahubs.net/download/project_code_donor

