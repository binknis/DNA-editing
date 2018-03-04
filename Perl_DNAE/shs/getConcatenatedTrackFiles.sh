
datadir=$1
org=$2 #
class=$3 #LTR
mm=$4 #GA
trackDirRegex=$5  #tracks*1e-0_5 #*
argPair=$6 #1e-0_5
filter=$7
outdir=$8 

mm=`echo -n $mm | tr '[a-z]' '[A-Z]'` 
mmS=${mm:0:1}
mmT=${mm:1:2}

# echo $datadir/$org/$class/results/Tracks/$trackDirRegex/$mm/$filter
# find $datadir/$org/$class/results/Tracks/$trackDirRegex/$mm/$filter
# exit;

outdir=$outdir/$filter/$class #added on 18/2/18
mkdir -p $outdir

for nuc in $mmS $mmT
do
	#nucListFreqInSeqPerPair
	find $datadir/$org/$class/results/Tracks/$trackDirRegex/$mm/$filter -name 'nucListFreqInSeqPerPair_'$nuc'*' -exec cat {} \; | sed -e 's/=/\t/g' | awk '{print $0 "\t" "'$nuc'"}' > $outdir/nucListFreqInSeqPerPair_"$nuc"_"$filter"_"$mm".txt
	
	#nucListFreqInConsPerPair
	find $datadir/$org/$class/results/Tracks/$trackDirRegex/$mm/$filter -name 'nucListFreqInConsPerPair_'$nuc'*' -exec cat {} \; | sed -e 's/=/\t/g' | awk '{print $0 "\t" "'$nuc'"}' > $outdir/nucListFreqInConsPerPair_"$nuc"_"$filter"_"$mm".txt
	
	# nucListFreq
	find $datadir/$org/$class/results/Tracks/$trackDirRegex/$mm/$filter -name 'nucListFreqInCons_'$nuc'*' -exec cat {} \; | sed -e 's/=/\t/g' | awk '{print $0 "\t" "'$nuc'"}' > $outdir/nucListFreqInCons_"$nuc"_"$filter"_"$mm".txt

	# motifPerSeq (1e-2) - 
	find $datadir/$org/$class/results/Tracks/$trackDirRegex/$mm/$filter -name 'motifPerSeq_'$nuc'_*1e-2*' -exec grep -H "" {} \; | sed -e 's/.*\/motifPerSeq_'$nuc'_//' -e 's/_1e-0_5.txt:/\t/' -e 's/_/\t/' -e 's/_/\t/' -e 's/_/\t/' -e 's/=/\t/g' | awk 'BEGIN{OFS="\t"} {print $5,$2,$3,$4,$9,$6,$10}' > $outdir/motifPerSeq_"$nuc"_1e-2_"$filter"_"$mm".txt

	# motifPerSeq (1e-2) - 
	find $datadir/$org/$class/results/Tracks/$trackDirRegex/$mm/$filter -name 'motifPerSeq_'$nuc'_*1e-3*' -exec grep -H "" {} \; | sed -e 's/.*\/motifPerSeq_'$nuc'_//' -e 's/_1e-0_5.txt:/\t/' -e 's/_/\t/' -e 's/_/\t/' -e 's/_/\t/' -e 's/=/\t/g' | awk 'BEGIN{OFS="\t"} {print $5,$2,$3,$4,$9,$6,$10}' > $outdir/motifPerSeq_"$nuc"_1e-3_"$filter"_"$mm".txt

	# nucComposition
	find $datadir/$org/$class/results/Tracks/$trackDirRegex/$mm/$filter -name 'nucComp_'$nuc'*' -exec grep -H "" {} \; | sed -e 's/.*\/nucComp_'$nuc'_//' -e 's/.*\/_//' -e 's/_1e-0_5.txt:/\t/' -e 's/_/\t/' -e 's/_/\t/' -e 's/=/\t/g' | awk '$4 != "ALL"' > $outdir/nucComp_"$nuc"_"$filter"_"$mm".txt

	# nucCompositionBordered
	find $datadir/$org/$class/results/Tracks/$trackDirRegex/$mm/$filter -name 'nucComp_bordered_'$nuc'*' -exec grep -H "" {} \; | sed -e 's/.*\/nucComp_bordered_'$nuc'_//' -e 's/.*\/_//' -e 's/_1e-0_5.txt:/\t/' -e 's/_/\t/' -e 's/_/\t/' -e 's/=/\t/g' | awk '$4 != "ALL"' > $outdir/nucComp_bordered_"$nuc"_"$filter"_"$mm".txt

	# sequence length 
	find $datadir/$org/$class/results/Tracks/$trackDirRegex/$mm/$filter -name 'siteList_'$nuc'*' -exec grep -H "" {} \; | sed -e 's/.*\/siteList_'$nuc'_//' -e 's/.*\/_//' -e 's/_1e-0_5.txt:/\t/' -e 's/_/\t/' -e 's/_/\t/' -e 's/=/\t/g' | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$5}' | sed -e 's/\S*:\([0-9]*\)-\([0-9]*\)[+-]$/\1\t\2/' | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$7-$6}' > $outdir/seqLen_"$nuc"_"$filter"_"$mm".txt

	# assemblyAndCoords
	find $datadir/$org/$class/results/Tracks/$trackDirRegex/$mm/$filter -name 'siteList_'$nuc'*' -exec cat {} \; | cut -f 1 | sed 's/=/\t/g' | cut -f 1-2 > $outdir/assemblyAndCoords_"$nuc"_"$filter"_"$mm".txt

	# First HSP in alignment of seqs to consensus (best in LTR class):
	find $datadir/$org/$class/results/Tracks/$trackDirRegex/$mm/$filter/SubfamFiles -name 'blastToClass_'$nuc'*' -exec cat {} \; |  awk 'a != $1 ; {a=$1}' > $outdir/blastToClass_"$nuc"_"$filter"_allOrgs_"$mm"_1stHSP.txt &

	#Num edited per
	find $datadir/$org/$class/results/Tracks/$trackDirRegex/$mm/$filter -name 'siteList_'$nuc'*' -exec cat {} \; | awk '{print $1 "\t" NF-1}' | sed 's/=/\t/g' | awk 'BEGIN{OFS="\t"}{print $1, $3, $4, $5, $2, $6}' > $outdir/numSitesPerSeq_"$nuc"_"$filter"_"$mm".txt

	#Freq (normalized freq)
	find $datadir/$org/$class/results/Tracks/$trackDirRegex/$mm/$filter -name 'freq_'$nuc'*'  -exec grep -H "" {} \; | sed -e 's/.*\///' -e 's/:/\t/' | sed 's/_/\t/g' | cut -f 3-5,8-12 > $outdir/freq_"$nuc"_"$filter"_"$mm".txt

	#Freq (normalized freq) 3range
	find $datadir/$org/$class/results/Tracks/$trackDirRegex/$mm/$filter -name 'freq_'$nuc'*'  -exec grep -H "" {} \; | sed -e 's/.*\///' -e 's/:/\t/' | awk '$2>-4 && $2<4' | sed 's/_/\t/g' | cut -f 3-5,8-12 > $outdir/freq_"$nuc"_"$filter"_3range.txt

	#Nuc-Frequencies per sequence
	find $datadir/$org/$class/results/Tracks/$trackDirRegex/$mm/$filter -name 'freqPerSeq_'$nuc'_*' -exec grep -H "" {} \; | sed -e 's/.*\/freqPerSeq_'$nuc'_//' -e 's/_/\t/' -e 's/_/\t/' -e 's/_.*\.txt:/\t/' | grep -vP '\tid\tpos\t' > $outdir/freqPerSeq_"$nuc"_"$filter"_"$mm".txt

	#Mapped consensi
	find $datadir/$org/$class/results/Tracks/*/$mm/$filter/SubfamFiles -name 'mapped_'$nuc'.txt'  -exec cat {} \; > $outdir/mapped_"$nuc"_"$filter"_"$mm".txt

	#posCons
	find $datadir/$org/$class/results/Tracks/*/$mm/$filter -name 'posCons_'$nuc'*' -exec grep -H \"\" {} \; | sed -e 's/.*\/posCons_'$nuc'_//' -e 's/.*\/_//' -e 's/_1e-0_5.txt:/\t/' -e 's/_/\t/' -e 's/_/\t/' -e 's/=/\t/g' > $outdir/posCons_"$nuc"_"$filter"_"$mm".txt	

done


#Non per nuc
# outDeg -
find $datadir/$org/$class/results/Tracks/$trackDirRegex/$mm/$filter -name 'outDeg*' -exec grep -H "" {} \; | sed -e 's/.*\/outDeg_//' -e 's/_1e-0_5.txt:/\t/' -e 's/_/\t/' -e 's/_/\t/' -e 's/=/\t/' | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$5,$4,$6}' > $outdir/outDeg_"$filter"_"$mm".txt
# inDeg - 
find $datadir/$org/$class/results/Tracks/$trackDirRegex/$mm/$filter -name 'inDeg*' -exec grep -H "" {} \; | sed -e 's/.*\/inDeg_//' -e 's/_1e-0_5.txt:/\t/' -e 's/_/\t/' -e 's/_/\t/' -e 's/=/\t/' | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$5,$4,$6}' > $outdir/inDeg_"$filter"_"$mm".txt
# Pairs (graph2 file)
find $datadir/$org/$class/results/Tracks/$trackDirRegex/$mm/$filter -name 'graph2*' -exec cat {} \; > $outdir/graph2_"$filter"_allOrgs_"$mm".txt &
