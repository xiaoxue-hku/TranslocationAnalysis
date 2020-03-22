# $1 cancer type

# $2 SV type

count()
{
	intersectBed -a $1_$2.bed -b ../src/hg19.knownCanonical+2kb.autosome.bed.txt -wo | sortBed -i - | uniq - | \
	cut -f12 - | sort - | uniq -c | awk -F ' ' '{print $2 "\t" $1}' - | \
	awk -F '\t' 'NR==FNR{id[$1]=$2;next} {if($2 in id) print $2 "\t" id[$2]; else print $2 "\t" "0"}' - ../src/convertID.tsv | sort - | uniq - > $1."$2"num.tsv
}

CountFPKM()
{
	awk -F '\t' 'NR !=1 {print (NF-1)*0.8}' ../src/TCGA-$1.htseq_fpkm-uq.tsv > $1.countCutoff

	awk -F '\t' 'NR !=1 {count=0; for(N=2; N<=NF; N++) {if($N==0) count++}; print count}' ../src/TCGA-$1.htseq_fpkm-uq.tsv > $1.countZero

	awk -F '\t' 'NR != 1 {T=0; for(N=2; N<=NF; N++) T+=$N; print T/=(NF-1)}' ../src/TCGA-$1.htseq_fpkm-uq.tsv > $1.countAverage

	paste ../src/convertID.tsv $1.countAverage $1.countCutoff $1.countZero | awk -F '\t' '{if($3==0 || $5 > $4) print $1 "\t" $2 "\t" "0"; else print $1 "\t" $2 "\t" $3}' - | sort -k2 - > $1.fpkm.tsv

	rm $1.countAverage $1.countCutoff $1.countZero
}

AnnotateSVFPKM()
{
	countSample=$(cut -f7 $1_$2.bed | sort - | uniq - | wc -l)	#7th column is donor id; check sample number by 'echo $countSample'

	awk -F '\t' 'NR==FNR{id[$2]=$2 "\t" $3;next} $1 in id {print id[$1] "\t" $2}' $1.fpkm.tsv ../src/hg19_HEPG2.rep.tsv | \
	awk -F '\t' 'NR==FNR{sv[$1]=$0;next} $1 in sv {print sv[$1] "\t" $2 "\t" $3}' $1."$2"num.tsv - | sort -k1 - | uniq - | \
	awk -F '\t' -v n=$countSample '{print $0 "\t" n}' - > $1."$2"Num.fpkm.repTime.sampleNum.tsv

	awk -F '\t' 'NR==FNR {a[$4]=$0; next} $1 in a {print a[$1] "\t" $0}' ../src/hg19.knownCanonical+2kb.autosome.bed.txt $1."$2"Num.fpkm.repTime.sampleNum.tsv | sortBed -i - > $1.annotated."$2"Num.fpkm.repTime.sampleNum.bed

	rm $1.fpkm.tsv $1."$2"num.tsv $1."$2"Num.fpkm.repTime.sampleNum.tsv
}


ComputeDistribution()
{
	countSample=$(cut -f7 $1_$2.bed | sort - | uniq - | wc -l)
	
	S=$(awk -F '\t' -v n=$countSample '{if($9==0) {SV+=$8; len+=($3-$2);}} END {print 1000000*SV/(len*n)}' "$1".annotated."$2"Num.fpkm.repTime.sampleNum.bed)

	NS=$(awk -F '\t' -v n=$countSample '{if($9>=1) {SV+=$8; len+=($3-$2);}} END {print 1000000*SV/(len*n)}' "$1".annotated."$2"Num.fpkm.repTime.sampleNum.bed)

	Innum=$(intersectBed -a "$1"_"$2".bed -b ../src/hg19.intergenic-2kb.autosome.bed.txt -u | wc -l)

	Inlen=$(awk -F '\t' '{len += ($3-$2)} END {print len}' ../src/hg19.intergenic-2kb.autosome.bed.txt)
	
	Intergenic=$(awk -v SV=$Innum -v len=$Inlen -v n=$countSample 'BEGIN {print 1000000*SV/(len*n)}')
	
	awk -v NS=$NS -v S=$S -v Intergenic=$Intergenic -v cancer=$1 'BEGIN {print cancer "\t" NS "\t" S "\t" Intergenic}' > $1.count$2.NSSIntergenic.tsv

	Rscript plotDistribution.r $1.count$2.NSSIntergenic.tsv
}

Compute5fpkmbs()
{
	awk -F '\t' '{if($9==0) print}' $1.annotated."$2"Num.fpkm.repTime.sampleNum.bed | sortBed -i - | uniq - > "$1"_"$2"bin04
	
	awk -F '\t' '{if($9!=0) print}' $1.annotated."$2"Num.fpkm.repTime.sampleNum.bed | sortBed -i - | uniq - | sort -k9nr - > $1_$2.Nzero
	
	N=$(wc -l < "$1"_"$2".Nzero)
	
	split -l $(expr $N / 4 + 1) -d $1_$2.Nzero "$1"_"$2"bin		#output "$1"_"$2"_"$3"bin0[0-3]
	
	for i in $(seq 0 4); do awk -F '\t' '{fpkm+=$9} END {print fpkm/NR}' "$1"_"$2"bin0"${i}" > Avgfpkm${i}; done
	
	for i in $(seq 0 4); do awk -F '\t' '{SV+=$8; len+=($3-$2);} END {print 1000000*SV/($11*len)}' "$1"_"$2"bin0"${i}" > AvgSV${i}; done
	
	AvgSV=$(paste -d "\n" AvgSV*)
	
	Avgfpkm=$(paste -d "\n" Avgfpkm*)
	
	paste -d "\t" <(printf "%s" "$Avgfpkm") <(printf "%s" "$AvgSV") | sort -k1n - > $1_$2.fpkmbs.tsv		#avgFPKM ~ avgSV
	
	rm AvgSV* Avgfpkm* "$1"_"$2"bin* "$1"_"$2".Nzero
	
	Rscript plotBinnedResults.r $1_$2.fpkmbs.tsv	
}


Compute5repTimebs()
{

	sort -k10nr "$1".annotated."$2"Num.fpkm.repTime.sampleNum.bed > "$1"_"$2".tmp
	
	N=$(wc -l < "$1"_"$2".tmp)
	
	split -l $(expr $N / 5 + 1) -d "$1"_"$2".tmp "$1"_"$2"repbin
	
	for i in $(seq 0 4); do awk -F '\t' '{rep+=$10} END {print rep/NR}' "$1"_"$2"repbin0"${i}" > Avgrep${i}; done
	
	for i in $(seq 0 4); do awk -F '\t' '{SV+=$8; len+=($3-$2);} END {print 1000000*SV/($11*len)}' "$1"_"$2"repbin0"${i}" > AvgSV${i}; done
	
	AvgSV=$(paste -d "\n" AvgSV*)
	
	Avgrep=$(paste -d "\n" Avgrep*)
	
	paste -d "\t" <(printf "%s" "$Avgrep") <(printf "%s" "$AvgSV") > $1_$2.repbs.tsv		#avgReplicationTiming ~ avgSV 
	
	rm AvgSV* Avgrep* "$1"_"$2"repbin* "$1"_"$2".tmp
	
	Rscript plotBinnedResults.r $1_$2.repbs.tsv	
}

ComputeProfile()
{
	ComputeCoverage() 
	{

		coverageBed -a $3 -b $1_$2.bed -d | python ../src/makeCoverageProfile.py - $4 $5

		tail -n 1 "$4"_binned.txt | cut -f2- - > $4.txt				

		#Regions are divided into silent and non-silent subgroups

		awk -F '\t' 'NR==FNR{if($9==0) id[$4]=$4;next} $4 in id{print $0}' $1.annotated."$2"Num.fpkm.repTime.sampleNum.bed $3 | sortBed -i - | uniq - > "$4"_S.bed.txt

		awk -F '\t' 'NR==FNR{if($9>=1) id[$4]=$4;next} $4 in id{print $0}' $1.annotated."$2"Num.fpkm.repTime.sampleNum.bed $3 | sortBed -i - | uniq - > "$4"_NS.bed.txt

		coverageBed -a "$4"_S.bed.txt -b $1_$2.bed -d | python ../src/makeCoverageProfile.py - "$4"_S $5

		coverageBed -a "$4"_NS.bed.txt -b $1_$2.bed -d | python ../src/makeCoverageProfile.py - "$4"_NS $5

		tail -n 1 "$4"_S_binned.txt | cut -f2- - > "$4"_S.txt

		tail -n 1 "$4"_NS_binned.txt | cut -f2- - > "$4"_NS.txt
		
	}

	ComputeCoverage $1 $2 ../src/hg19.tss-2kb.bed.txt "$1"_"$2"_tss-2kb 20

	ComputeCoverage $1 $2 ../src/hg19.knownCanonical.autosome.bed.txt "$1"_"$2"_body 50

	ComputeCoverage $1 $2 ../src/hg19.tes+2kb.bed.txt "$1"_"$2"_tes+2kb 20

	paste -d "" "$1"_"$2"_tss-2kb.txt "$1"_"$2"_body.txt "$1"_"$2"_tes+2kb.txt > "$1"_"$2".all.profile			# SVNum / Mb

	paste -d "" "$1"_"$2"_tss-2kb_S.txt "$1"_"$2"_body_S.txt "$1"_"$2"_tes+2kb_S.txt > "$1"_"$2".S.profile

	paste -d "" "$1"_"$2"_tss-2kb_NS.txt "$1"_"$2"_body_NS.txt "$1"_"$2"_tes+2kb_NS.txt > "$1"_"$2".NS.profile

	rm *"$1"_"$2"*.txt
	
	Rscript plotProfile.r "$1"_"$2".all.profile
	
	Rscript plotProfile.r "$1"_"$2".NS.profile
	
	Rscript plotProfile.r "$1"_"$2".S.profile
}

TranslocationHR()
{
	awk -F '\t' 'NR==FNR{id[$1]=$4;next} $7 in id {print $0 "\t" id[$7]}' ../src/PCAWG_sub_signatures_in_samples_beta2.20170320.donor $1_$2.bed | \
	awk -F '\t' '{if($9==0) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"HRpos"}' - | sortBed -i - > "$1"HRpos_$2.bed
	
	awk -F '\t' 'NR==FNR{id[$1]=$4;next} $7 in id {print $0 "\t" id[$7]}' ../src/PCAWG_sub_signatures_in_samples_beta2.20170320.donor $1_$2.bed | \
	awk -F '\t' '{if($9!=0) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"HRpos"}' - | sortBed -i - > "$1"HRneg_$2.bed
	
	mv ../src/TCGA-$1.htseq_fpkm-uq.tsv ../src/TCGA-"$1"HRpos.htseq_fpkm-uq.tsv
	count "$1"HRpos $2
	CountFPKM "$1"HRpos $2
	AnnotateSVFPKM "$1"HRpos $2
	ComputeDistribution "$1"HRpos $2
	ComputeProfile "$1"HRpos $2
	
	mv ../src/TCGA-"$1"HRpos.htseq_fpkm-uq.tsv ../src/TCGA-"$1"HRneg.htseq_fpkm-uq.tsv
	count "$1"HRneg $2
	CountFPKM "$1"HRneg $2
	AnnotateSVFPKM "$1"HRneg $2
	ComputeDistribution "$1"HRneg $2	
	ComputeProfile "$1"HRneg $2
	
	mv ../src/TCGA-"$1"HRneg.htseq_fpkm-uq.tsv ../src/TCGA-$1.htseq_fpkm-uq.tsv
}

ComputeRegression()
{
	awk -F '\t' '{print $4 "\t" 1000000*$8/(($3-$2)*$11) "\t" $9 "\t" $10}' $1.annotated."$2"Num.fpkm.repTime.sampleNum.bed | sort -k1 - | \
	awk -F '\t' 'BEGIN{print "ID\ttranslocation/Mb/sample\tgene_expression\treplication_timing"}; {print}' - > ../regression/all/"$1"_"$2".all_regression.tsv

	sed 1d ../regression/all/"$1"_"$2".all_regression.tsv | awk -F '\t' '{if($3==0) print}' - | \
	awk -F '\t' 'BEGIN{print "ID\ttranslocation/Mb/sample\tgene_expression\treplication_timing"}; {print}' - > ../regression/silent/"$1"_"$2".S_regression.tsv
	
	sed 1d ../regression/all/"$1"_"$2".all_regression.tsv | awk -F '\t' '{if($3>=1) print}' - | \
	awk -F '\t' 'BEGIN{print "ID\ttranslocation/Mb/sample\tgene_expression\treplication_timing"}; {print}' - > ../regression/nonsilent/"$1"_"$2".NS_regression.tsv
	
	Rscript regression_on_fpkm_repTime.r ../regression/all/"$1"_"$2".all_regression.tsv
	
	Rscript regression_on_fpkm_repTime.r ../regression/silent/"$1"_"$2".S_regression.tsv
	
	Rscript regression_on_fpkm_repTime.r ../regression/nonsilent/"$1"_"$2".NS_regression.tsv
}

IntergenicRegression()
{
	cp $1_$2.bed ../regression && cd ../regression 

	mv ../src ../srctmp
	mv ../regressionSrc ../src
	mv ../src/hg19_numbered.intergenic-2kb.autosome.bed.txt ../src/hg19.knownCanonical+2kb.autosome.bed.txt 
	mv ../src/convertID_intergenic.tsv ../src/convertID.tsv 
	mv ../src/hg19_HEPG2_intergenic.rep.tsv ../src/hg19_HEPG2.rep.tsv
	
	count $1 $2
	awk -F '\t' '{print $0 "\t" "NA"}' ../src/convertID.tsv > $1.fpkm.tsv
	AnnotateSVFPKM $1 $2
	
	awk -F '\t' '{print $4 "\t" 1000000*$8/(($3-$2)*$11) "\t" $9 "\t" $10}' $1.annotated."$2"Num.fpkm.repTime.sampleNum.bed | sort -k1 - | \
	awk -F '\t' 'BEGIN{print "ID\ttranslocation/Mb/sample\tgene_expression\treplication_timing"}; {print}' - > ./intergenic/"$1"_"$2".intergenic_regression.tsv
	
	mv ../src/hg19.knownCanonical+2kb.autosome.bed.txt ../src/hg19_numbered.intergenic-2kb.autosome.bed.txt 
	mv ../src/convertID.tsv ../src/convertID_intergenic.tsv 
	mv ../src/hg19_HEPG2.rep.tsv ../src/hg19_HEPG2_intergenic.rep.tsv 
	mv ../src ../regressionSrc
	mv ../srctmp ../src
	
	rm *.bed && cd ../$2
	
	Rscript regression_on_fpkm_repTime.r ../regression/intergenic/"$1"_"$2".intergenic_regression.tsv
}

HRregression()
{
	awk -F '\t' '{print $0 "\t" "+"}' "$1"HRpos.annotated."$2"Num.fpkm.repTime.sampleNum.bed > $1_$2.HRpos
	
	awk -F '\t' '{print $0 "\t" "-"}' "$1"HRneg.annotated."$2"Num.fpkm.repTime.sampleNum.bed > $1_$2.HRneg
	
	cat $1_$2.HRpos $1_$2.HRneg | sortBed -i - | awk -F '\t' '{print $4 "\t" 1000000*$8/(($3-$2)*$11) "\t" $9 "\t" $10 "\t" $12}' - | sort -k1 - | \
	awk -F '\t' 'BEGIN{print "ID\ttranslocation/Mb/sample\tgene_expression\treplication_timing\tHR"}; {print}' - > ../regression/HRall/"$1"_"$2".HRall_regression.tsv
	
	sed 1d ../regression/HRall/"$1"_"$2".HRall_regression.tsv | awk -F '\t' '{if($3>=1) print}' - | \
	awk -F '\t' 'BEGIN{print "ID\ttranslocation/Mb/sample\tgene_expression\treplication_timing\tHR"}; {print}' - > ../regression/HRnonsilent/"$1"_"$2".HRNS_regression.tsv

	rm $1_$2.HRpos $1_$2.HRneg
	
	Rscript regression_on_fpkm_repTime_HR.r ../regression/HRall/"$1"_"$2".HRall_regression.tsv
	
	Rscript regression_on_fpkm_repTime_HR.r ../regression/HRnonsilent/"$1"_"$2".HRNS_regression.tsv
}

oneMbRepTimeregression()
{
	cp $1_$2.bed ../regression && cd ../regression 

	mv ../src ../srctmp
	mv ../regressionSrc ../src
	mv ../src/hg19_numbered.1Mb.autosome.bed.txt ../src/hg19.knownCanonical+2kb.autosome.bed.txt 
	mv ../src/convertID_1Mb.tsv ../src/convertID.tsv 
	mv ../src/hg19_HEPG2_1Mb.rep.tsv ../src/hg19_HEPG2.rep.tsv
	
	count $1 $2
	awk -F '\t' '{print $0 "\t" "NA"}' ../src/convertID.tsv > $1.fpkm.tsv
	AnnotateSVFPKM $1 $2
	Compute5repTimebs $1 $2 && mv $1_$2.repbs.tsv ./rep1Mb/$1_$2_.1Mbrepbs.tsv
	
	awk -F '\t' '{print $4 "\t" 1000000*$8/(($3-$2)*$11) "\t" $9 "\t" $10}' $1.annotated."$2"Num.fpkm.repTime.sampleNum.bed | sort -k1 - | \
	awk -F '\t' 'BEGIN{print "ID\ttranslocation/Mb/sample\tgene_expression\treplication_timing"}; {print}' - > ./rep1Mb/"$1"_"$2".1Mb_regression.tsv
	
	mv ../src/hg19.knownCanonical+2kb.autosome.bed.txt ../src/hg19_numbered.1Mb.autosome.bed.txt 
	mv ../src/convertID.tsv ../src/convertID_1Mb.tsv
	mv ../src/hg19_HEPG2.rep.tsv ../src/hg19_HEPG2_1Mb.rep.tsv
	mv ../src ../regressionSrc
	mv ../srctmp ../src
	
	rm *.bed && cd ../$2
	
	Rscript regression_on_repTime_for_1Mb_bin ../regression/rep1Mb/"$1"_"$2".1Mb_regression.tsv
}

