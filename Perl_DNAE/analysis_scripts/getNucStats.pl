#Function: Create all the nuc-stat files for a seqFasta and siteList file pair
#Input: a directory that contains a seqFasta and siteList file for a specific organism 
#	Notes: the files must have same name and only differ in seqFasta -> siteList difference and extension of .fa -> .txt (examples: )

use strict;
use lib $ENV{"HOME"} ."/Perl_DNAE";
use getAll; 
use lib $ENV{"HOME"} ."/Perl_DNAE/analysis_scripts"; 
use analysisSubs; 
use File::Path qw(mkpath);
my $home = $ENV{"HOME"}; 


(my $dir, my $GA, my $siteListFile) = @ARGV; 
my $suffix = $siteListFile;
$suffix =~ s/.*siteList_//; 
$suffix =~ s/\.txt$//; 

$dir .= "/"; 
$GA = "G" unless $GA; #default "G"
my $ga = lc $GA; 
my $trim = 0; #if to trim polyA tail using trimest
my $border = 0; #if to calc nuc freq only within border of first and last editing sites #***
my $nucleotides = 0; #"cgt"; 

# my $seqFile = $dir . "seqFasta_".$GA."_Human_SINE_Alu_1e-0_4.fa";
my $seqFile = $dir . "seqFasta_". $suffix .".fa";
if ($trim){
	my $trimmedSeqFile = $seqFile; 
	$trimmedSeqFile =~ s/seqFasta/seqFastaTrimmed/; 
	my $mismatches = 1; 
	analysisSubs::trimPolyA($seqFile, $trimmedSeqFile, 0, $mismatches);
	$seqFile = $trimmedSeqFile; 
}

$siteListFile = $dir . "siteList_".$suffix.".txt"; 
# my $range = 2; 
# my $freqFile = $dir . "logo".($range==2 ? "" : $range)."_".$suffix."_freq.txt";
my $range = 7; 
my $freqFile = $dir . "rawFreq_".$suffix.".txt";
my $freqPerSeqFile = $dir . "freqPerSeq_".$suffix.".txt"; 
my $no_CpG = 0; 
analysisSubs::getNucFrequencyPerPosAllSeqs($seqFile, $siteListFile, $range, 1, 0, $freqFile, $no_CpG, $freqPerSeqFile); #(my $seqFile, my $siteListFile, my $range, my $normalize, my $reverse_editing, my $outfile, my $no_CpG, my $freqPerSeqFile)

my $hash_ref; 
if ($border){
	$hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $seqFile, $siteListFile, $nucleotides); #norm within borders!!!
}
else {
	$hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $seqFile, 0, $nucleotides);
}
 my $freqNormedFile = $siteListFile; 
 $freqNormedFile =~ s/.*siteList/freq/; 
 $freqNormedFile = $dir .($border ? "bordered_" : "").($nucleotides ? $nucleotides ."_" : ""). $freqNormedFile; 
analysisSubs::printFreqHash($hash_ref, $freqNormedFile);

#print background frequencies used for above normalization
my $retfmt = 1;
my $bgFreq_file = $freqNormedFile; 
$bgFreq_file =~ s/freq/bgFreq/;
my $siteListParam = ($border ? $siteListFile : 0); 
my $bg_hash_ref = analysisSubs::getNucFreqFromFasta($seqFile, $range, $retfmt, $siteListParam);
analysisSubs::printFreqHash($bg_hash_ref->{$ga}, $bgFreq_file);

#print background for chi-squared test (background frequencies - edited position frequencies) 
#Explanation: needed because contingency tables need disjoint groups and the previous "background" contains the edited positions too. This doesn't. 
my $normalize = 1; my $seqFile2 = 0; my $reverse_editing=0; 
my $bgneFreq_file = $freqNormedFile; 
$bgneFreq_file =~ s/freq/bgNoEditedFreq/;
my $bgne_hash_ref = analysisSubs::getBackgroundNoEditedFreq($GA, $seqFile, $siteListFile, $range, $seqFile2, $normalize, $bgneFreq_file, $reverse_editing); #function also prints so don't need to do so here

#get Triplets (enrichment)
my $acgt = $ga; my $normalize = $seqFile; my $suppressPrint = 1; my $revcom = 0; my $range=1; $reverse_editing=0;  
# my $acgt = $ga; my $normalize = $seqFile; my $suppressPrint = 0; my $revcom = 0; my $range=1; my $reverse_editing=0; my $border=1; #norm within borders!!!
analysisSubs::getTripletsAllSeqs($seqFile, $acgt, $siteListFile, $range, $normalize, $reverse_editing, $suppressPrint, $revcom, $border);
#get Triplets (frequencies)
$normalize = 1;
analysisSubs::getTripletsAllSeqs($seqFile, $acgt, $siteListFile, $range, $normalize, $reverse_editing, $suppressPrint, $revcom, $border);

#get Triplets (enrichment)
my $acgt = $ga; my $normalize = $seqFile; my $suppressPrint = 1; my $revcom = 0; my $range=1; $reverse_editing=0;  
# my $acgt = $ga; my $normalize = $seqFile; my $suppressPrint = 0; my $revcom = 0; my $range=1; my $reverse_editing=0; my $border=1; #norm within borders!!!
analysisSubs::getTripletsAllSeqs($seqFile, $acgt, $siteListFile, $range, $normalize, $reverse_editing, $suppressPrint, $revcom, $border);
#get Triplets (frequencies)
$normalize = 1;
analysisSubs::getTripletsAllSeqs($seqFile, $acgt, $siteListFile, $range, $normalize, $reverse_editing, $suppressPrint, $revcom, $border);

#get quintuplets (frequencies and fractions)
my $acgt = $ga; my $normalize = $seqFile; my $suppressPrint = 1; my $revcom = 0; my $range=2; $reverse_editing=0;  
$normalize = 0;
analysisSubs::getTripletsAllSeqs($seqFile, $acgt, $siteListFile, $range, $normalize, $reverse_editing, $suppressPrint, $revcom, $border);
$normalize = 1;
analysisSubs::getTripletsAllSeqs($seqFile, $acgt, $siteListFile, $range, $normalize, $reverse_editing, $suppressPrint, $revcom, $border);
#print quintuplet BG (fraction) to file
my $tripletBGfile = $seqFile; 
$normalize = 1; 
my $triBGprefix = $normalize ? "neighbors2BGfrac" : "neighbors2BG"; 
$tripletBGfile =~ s/seqFasta/$triBGprefix/; 
$tripletBGfile =~ s/\.fa$/.txt/;
my $siteListArg = 0; 
analysisSubs::printTripletBG($seqFile, $acgt, $siteListArg, $range, $normalize, $tripletBGfile, $revcom);
#print quintuplet BG (count) to file
my $tripletBGfile = $seqFile; 
$normalize = 0; 
my $triBGprefix = $normalize ? "neighbors2BGfrac" : "neighbors2BG"; 
$tripletBGfile =~ s/seqFasta/$triBGprefix/; 
$tripletBGfile =~ s/\.fa$/.txt/;
my $siteListArg = 0; 
analysisSubs::printTripletBG($seqFile, $acgt, $siteListArg, $range, $normalize, $tripletBGfile, $revcom);

#create nucListFreq file
# my $nucListFile = $dir . "nucList_".$suffix.".txt"; 
# analysisSubs::nucListToFreq($nucListFile);

#create nucList + nucListFreq per pairs
# (my $suffix2) = $suffix =~ /[ACTG]_(\S+)/;
# my $clustersFile = $dir . "clusters_".$suffix2.".tab"; 
# my $SorT = ($GA =~ /[AT]/ ? 1 : 0); 
# analysisSubs::nucFreqPerPairs($clustersFile, $SorT, $GA);

#print nuc composition for all seqs and per seq
# my $ncFile = $seqFile; 
# $ncFile =~ s/seqFasta/nucComp/; 
# $ncFile =~ s/\.fa$/.txt/; 
# my $ncPerSeqFile = $seqFile; 
# $ncPerSeqFile =~ s/seqFasta/nucCompPerSeq/; 
# $ncPerSeqFile =~ s/\.fa$/.txt/;
###Note: Fields for printAllSeqNucComposition (my $seqFile, my $retfmt, my $siteListFile, my $outFile, my $outFilePerSeq)
# analysisSubs::printAllSeqNucComposition($seqFile,1,0,$ncFile, $ncPerSeqFile); #no bordering 
# $ncFile =~ s/nucComp/nucComp_bordered/; 
# $ncPerSeqFile =~ s/nucCompPerSeq/nucCompPerSeq_bordered/; 
# analysisSubs::printAllSeqNucComposition($seqFile, 1, $siteListFile, $ncFile, $ncPerSeqFile); #composition between terminal editing sites