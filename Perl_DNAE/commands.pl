use strict; 
use getAll; 
use lib "analysis_scripts"; 
use analysisSubs;


#Create file of nucComposition for DNA consensus sequences
my $seqFile = "/home/alu/binknis/binknis_data/RepBase/ConsForMapping_12_1_14/DNA.fa";
my $retfmt = 1; 
my $siteListFile = 0; 
my $outFile = "/home/alu/binknis/binknis_data/RepBase/nucComposition_DNA_fromConsForMapping_12_1_14.txt"; 
analysisSubs::printAllSeqNucComposition($seqFile, $retfmt, $siteListFile, $outFile); 

exit; 

#Create file of nucComposition for LTR consensus sequences
my $seqFile = "/home/alu/binknis/binknis_data/RepBase/ConsForMapping_12_1_14/LTR.fa";
my $retfmt = 1; 
my $siteListFile = 0; 
my $outFile = "/home/alu/binknis/binknis_data/RepBase/nucComposition_LTR_fromConsForMapping_12_1_14.txt"; 
analysisSubs::printAllSeqNucComposition($seqFile, $retfmt, $siteListFile, $outFile); 

exit; 

#For Natasa and Nika: align Lizard L1 to consensus to see editing sites
my $siteListFile = "/home/alu/binknis/Data/Lizard/LINE/results/Tracks/tracks_Lizard_LINE_L1_1e-8_8/siteList_A_Lizard_LINE_L1_1e-8_8.txt"; 
analysisSubs::getEditedPositionsInCons($siteListFile); 

exit; 

#check getTerminalsFromSiteList
my $siteListFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-0_4_iff4sites_oldANewGToAllPrimates_0.75uniqueAMap2A_bestSourcesIncomplete/siteList_A_Human_SINE_Alu_1e-0_4.txt";
my $tSites  = analysisSubs::getTerminalsFromSiteList($siteListFile);
foreach my $k(keys %$tSites){
	print $k . "\t".${$tSites->{$k}}[0] . "\t". ${$tSites->{$k}}[1]. "\n"; 
}
exit; 



#blast A and G vs consensus as 3-base
my $dir = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-0_5/";
my $qFa = $dir ."seqFasta_A_Human_SINE_Alu_1e-0_5.fa"; 
my $sFa = "/home/alu/binknis/binknis_data/RepBase/ConsForMapping_18_8_13/SINE/Alu.fa"; 
my $blastOut = $dir . "seqA_blast3baseVsAllCons.txt";
my $transform = "ga"; 
my $outfmt = 6; 
my $numThreads = 11; 
my $maxTargets = 5; 
analysisSubs::blast3base($qFa, $sFa, $blastOut,  $transform,  $outfmt,  $numThreads,  $maxTargets);

my $qFa2 = $qFa; 
$qFa2 =~ s/_A_/_G_/; 
$blastOut = $dir . "seqG_blast3baseVsAllCons.txt";
analysisSubs::blast3base($qFa2, $sFa, $blastOut, $transform,  $outfmt,  $numThreads,  $maxTargets); 

exit; 

#test getPairSites
my $clustFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped_bestPairs/clusters_Human_SINE_Alu_1e-0_5.tab"; 
my $ps = analysisSubs::getPairSites($clustFile); 
foreach my $k (keys %$ps){
	print $k ,"\t", join(',', @{$ps->{$k}}) ,"\n"; 
}
exit;

#trimmed file  - freq and trios
my $dir = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-0_5_5mapToA_allOfMapped_bestPairs_oldNew/";
my $seqFile = $dir . "seqFasta_A_Human_SINE_Alu_1e-0_5.fa";
my $trimmedSeqFile = $dir . "seqFastaTrimmed_A_Human_SINE_Alu_1e-0_5.fa";
my $mismatches = 1; 
analysisSubs::trimPolyA($seqFile, $trimmedSeqFile, 0, $mismatches); 

my $siteListFile = $dir . "siteList_A_Human_SINE_Alu_1e-0_5.txt"; 
my $freqFile = $dir . "logo_A_Human_SINE_Alu_1e-0_5_freq.Trimmed.txt"; 
my $range = 2; 
analysisSubs::getNucFrequencyPerPosAllSeqs($trimmedSeqFile, $siteListFile, $range, 1, 0, $freqFile); #(my $seqFile, my $siteListFile, my $range, my $normalize, my $reverse_editing, my $outfile)

# my $hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $trimmedSeqFile); 
my $hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $trimmedSeqFile, $siteListFile); #norm within borders!!!
analysisSubs::printFreqHash($hash_ref);

analysisSubs::countSeqSitesInTrimmed($trimmedSeqFile, $siteListFile);

# my $acgt = "a"; my $normalize = $seqFile; my $suppressPrint = 0; my $revcom = 0; my $range=1; my $reverse_editing=0; my $border=0; 
my $acgt = "a"; my $normalize = $seqFile; my $suppressPrint = 0; my $revcom = 0; my $range=1; my $reverse_editing=0; my $border=1; #norm within borders!!!
analysisSubs::getTripletsAllSeqs($seqFile, $acgt, $siteListFile, $range, $normalize, $reverse_editing, $suppressPrint, $revcom, $border);

exit; 

#not trimmed - freq and trios

my $dir = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-0_5_5mapToA_allOfMapped_bestPairs/";
my $seqFile = $dir . "seqFasta_A_Human_SINE_Alu_1e-0_5.fa";

my $siteListFile = $dir . "siteList_A_Human_SINE_Alu_1e-0_5.txt"; 
my $freqFile = $dir . "logo_A_Human_SINE_Alu_1e-0_5_freq.txt"; 
my $range = 2; 
analysisSubs::getNucFrequencyPerPosAllSeqs($seqFile, $siteListFile, $range, 1, 0, $freqFile); #(my $seqFile, my $siteListFile, my $range, my $normalize, my $reverse_editing, my $outfile)

# my $hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $seqFile); 
my $hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $seqFile, $siteListFile); #norm within borders!!!
analysisSubs::printFreqHash($hash_ref);

# my $acgt = "a"; my $normalize = $seqFile; my $suppressPrint = 0; my $revcom = 0; my $range=1; my $reverse_editing=0; my $border=0; 
my $acgt = "a"; my $normalize = $seqFile; my $suppressPrint = 0; my $revcom = 0; my $range=1; my $reverse_editing=0; my $border=1; #norm within borders!!!
analysisSubs::getTripletsAllSeqs($seqFile, $acgt, $siteListFile, $range, $normalize, $reverse_editing, $suppressPrint, $revcom, $border);

exit; 


#many analysis tracks_Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped_oldNew_bestPairs
my $dir = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped_oldNew_bestPairs/"; 
my $seqFile = $dir . "seqFasta_A_Human_SINE_Alu_1e-0_5.fa";
my $trimmedSeqFile = $dir . "seqFastaTrimmed_A_Human_SINE_Alu_1e-0_5.fa";
analysisSubs::trimPolyA($seqFile, $trimmedSeqFile); 

my $siteListFile = $dir . "siteList_A_Human_SINE_Alu_1e-0_5.txt"; 
my $freqFile = $dir . "logo_A_Human_SINE_Alu_1e-0_5_freq.Trimmed.txt"; 
my $range = 2; 
analysisSubs::getNucFrequencyPerPosAllSeqs($trimmedSeqFile, $siteListFile, $range, 1, 0, $freqFile); #(my $seqFile, my $siteListFile, my $range, my $normalize, my $reverse_editing, my $outfile)

#my $hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $trimmedSeqFile); 
my $hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $trimmedSeqFile, $siteListFile); #norm within borders!!!
analysisSubs::printFreqHash($hash_ref);

analysisSubs::countSeqSitesInTrimmed($trimmedSeqFile, $siteListFile);

#my $acgt = "a"; my $normalize = $seqFile; my $suppressPrint = 0; my $revcom = 0; my $range=1; my $reverse_editing=0; my $border=0; 
my $acgt = "a"; my $normalize = $seqFile; my $suppressPrint = 0; my $revcom = 0; my $range=1; my $reverse_editing=0; my $border=0; #norm within borders!!!
analysisSubs::getTripletsAllSeqs($seqFile, $acgt, $siteListFile, $range, $normalize, $reverse_editing, $suppressPrint, $revcom, $border);

exit; 

#best pairs 
my $dir = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped_oldNew/"; 
my $pairFile = $dir ."graph2_Human_SINE_Alu_1e-0_5.txt";
my $fastaS = $dir . "seqFasta_A_Human_SINE_Alu_1e-0_5.fa"; 
my $fastaT = $dir . "seqFasta_G_Human_SINE_Alu_1e-0_5.fa"; 
my $STreversed = 1; 
my $transform = "ga"; 
my $hash_ref = analysisSubs::getBestSources($pairFile, $fastaS, $fastaT, $STreversed, $transform); 
exit; 

#calc motif for Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped with polyA trimmed (trimmed regions discarded for editing sites and for normalization) - A seqs
my $seqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped/seqFasta_A_Human_SINE_Alu_1e-0_5.fa";
my $trimmedSeqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped/seqFastaTrimmed_A_Human_SINE_Alu_1e-0_5.fa";
analysisSubs::trimPolyA($seqFile, $trimmedSeqFile); 

my $siteListFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped/siteList_A_Human_SINE_Alu_1e-0_5.txt"; 
my $freqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped/logo_A_Human_SINE_Alu_1e-0_5_freq.Trimmed.txt"; 
my $range = 2; 
analysisSubs::getNucFrequencyPerPosAllSeqs($trimmedSeqFile, $siteListFile, $range, 1, 0, $freqFile); #(my $seqFile, my $siteListFile, my $range, my $normalize, my $reverse_editing, my $outfile)

my $hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $trimmedSeqFile); 
analysisSubs::printFreqHash($hash_ref);

analysisSubs::countSeqSitesInTrimmed($trimmedSeqFile, $siteListFile);

my $acgt = "a"; my $normalize = $seqFile; my $suppressPrint = 0; my $revcom = 0; my $range=1; my $reverse_editing=0; my $border=0; 
analysisSubs::getTripletsAllSeqs($seqFile, $acgt, $siteListFile, $range, $normalize, $reverse_editing, $suppressPrint, $revcom, $border);

exit; 

#Same as above for G-seqs: calc motif for Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped_oldNew with polyA trimmed (trimmed regions discarded for editing sites and for normalization)

my $seqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped/seqFasta_G_Human_SINE_Alu_1e-0_5.fa";
my $trimmedSeqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped/seqFastaTrimmed_G_Human_SINE_Alu_1e-0_5.fa";
analysisSubs::trimPolyA($seqFile, $trimmedSeqFile); 

my $siteListFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped/siteList_G_Human_SINE_Alu_1e-0_5.txt"; 
my $freqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped/logo_G_Human_SINE_Alu_1e-0_5_freq.Trimmed.txt"; 
my $range = 2; 
analysisSubs::getNucFrequencyPerPosAllSeqs($trimmedSeqFile, $siteListFile, $range, 1, 0, $freqFile); #(my $seqFile, my $siteListFile, my $range, my $normalize, my $reverse_editing, my $outfile)

my $hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $trimmedSeqFile); 
analysisSubs::printFreqHash($hash_ref);

analysisSubs::countSeqSitesInTrimmed($trimmedSeqFile, $siteListFile);

exit; 

#count number seqs + sites retained after trimming
my $siteListFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped_oldNew/siteList_A_Human_SINE_Alu_1e-0_5.txt"; 
my $trimmedSeqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped_oldNew/seqFastaTrimmed_A_Human_SINE_Alu_1e-0_5.fa";
analysisSubs::countSeqSitesInTrimmed($trimmedSeqFile, $siteListFile);
exit; 

#calc motif for Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped_oldNew with polyA trimmed (trimmed regions discarded for editing sites and for normalization)
my $seqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped_oldNew/seqFasta_A_Human_SINE_Alu_1e-0_5.fa";
my $trimmedSeqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped_oldNew/seqFastaTrimmed_A_Human_SINE_Alu_1e-0_5.fa";
analysisSubs::trimPolyA($seqFile, $trimmedSeqFile); 

my $siteListFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped_oldNew/siteList_A_Human_SINE_Alu_1e-0_5.txt"; 
my $freqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped_oldNew/logo_A_Human_SINE_Alu_1e-0_5_freq.Trimmed.txt"; 
my $range = 2; 
analysisSubs::getNucFrequencyPerPosAllSeqs($trimmedSeqFile, $siteListFile, $range, 1, 0, $freqFile); #(my $seqFile, my $siteListFile, my $range, my $normalize, my $reverse_editing, my $outfile)

my $hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $trimmedSeqFile); 
analysisSubs::printFreqHash($hash_ref);

exit; 
#calc enrichment scores for best pairs in 1e-0 5 (subset specified)
my $freqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped_oldNew/logo_A_Human_SINE_Alu_1e-0_5_freq.txt"; 
my $seqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped_oldNew/seqFasta_A_Human_SINE_Alu_1e-0_5.fa";
my $hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $seqFile); 
analysisSubs::printFreqHash($hash_ref);

exit; 

my $dir = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped/"; 
my $pairFile = $dir ."graph2_Human_SINE_Alu_1e-0_5.txt"; 
my $fastaS = $dir . "seqFasta_A_Human_SINE_Alu_1e-0_5.fa"; 
my $fastaT = $dir . "seqFasta_G_Human_SINE_Alu_1e-0_5.fa"; 
my $STreversed = 1; 
my $transform = "ga"; 
my $hash_ref = analysisSubs::getBestSources($pairFile, $fastaS, $fastaT, $STreversed, $transform); 
exit; 

#calc enrichment scores for best pairs in 1e-0 5 (subset specified)
my $freqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped/logo_A_Human_SINE_Alu_1e-0_5_freq.txt"; 
my $seqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped/seqFasta_A_Human_SINE_Alu_1e-0_5.fa";
my $hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $seqFile); 
analysisSubs::printFreqHash($hash_ref);

my $siteListFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-0_5_4mapToA_over0.8ofMapped/siteList_A_Human_SINE_Alu_1e-0_5.txt";
analysisSubs::getEditedPositionsInCons($siteListFile); 

exit; 
#get best A-sources for Alu 1e-0 5
my $dir = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-0_5/"; 
my $pairFile = $dir ."graph2_Human_SINE_Alu_1e-0_5.txt"; 
my $fastaS = $dir . "seqFasta_A_Human_SINE_Alu_1e-0_5.fa"; 
my $fastaT = $dir . "seqFasta_G_Human_SINE_Alu_1e-0_5.fa"; 
my $STreversed = 1; 
my $transform = "ga"; 
my $hash_ref = analysisSubs::getBestSources($pairFile, $fastaS, $fastaT, $STreversed, $transform); 
exit; 

#get count of each editing site mapped to consensus
# my $nucListFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-0_5/nucList_A_Human_SINE_Alu_1e-0_5.txt";
my $nucListFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-0_5/nucList_A_Human_SINE_Alu_1e-0_5.txt";
analysisSubs::nucListToFreq($nucListFile); 
exit; 

#get histogram of editing sites mapped to consensus
# my $siteListFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-0_5/siteList_G_Human_SINE_Alu_1e-0_5.txt";
# my $siteListFile = "/home/alu/binknis/Data/Human/LTR/results/Tracks/tracks_Human_LTR_ERVK_1e-7_7/siteList_G_Human_LTR_ERVK_1e-7_7.txt";
# my $siteListFile = "/home/alu/binknis/Data/Human/LINE/results/Tracks/tracks_Human_LINE_L1_1e-5_5/siteList_G_Human_LINE_L1_1e-5_5.txt"; 
# my $siteListFile = "/home/alu/binknis/Data/Chimp/SINE/results/Tracks/tracks_Chimp_SINE_Alu_1e-0_5/siteList_G_Chimp_SINE_Alu_1e-0_5.txt"; 
my $siteListFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-0_5/siteList_A_Human_SINE_Alu_1e-0_5.txt";
analysisSubs::getEditedPositionsInCons($siteListFile); 

exit; 
#get triple freqs for 1e-8 8 bestpairs
my $seqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-8_8_bestASources/seqFasta_A_Human_SINE_Alu_1e-8_8.fa";
my $acgt = "a"; 
my $siteListFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-8_8_bestASources/siteList_A_Human_SINE_Alu_1e-8_8.txt";
my $normalize = $seqFile; 
my $suppressPrint = 0;
my $revcom = 0;
my $range=1; 
my $reverse_editing=0; 
my $border=0; 
my $hash_ref = analysisSubs::getTripletsAllSeqs($seqFile, $acgt, $siteListFile, $range, $normalize, $reverse_editing, $suppressPrint, $revcom, $border);
exit; 

#normed nucFreq by whole A-seqs for 1e-0 6 
my $dir = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-0_6/"; 
my $freqFile = $dir . "logo_A_Human_SINE_Alu_1e-0_6_freq.txt"; 
my $seqFile = $dir ."seqFasta_A_Human_SINE_Alu_1e-0_6.fa"; 
my $hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $seqFile);
analysisSubs::printFreqHash($hash_ref);
exit; 

#get triple freqs for 1e-8 8 bestpairs
my $seqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-8_8_bestASources/seqFasta_A_Human_SINE_Alu_1e-8_8.fa";
my $acgt = "a"; 
my $siteListFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-8_8_bestASources/siteList_A_Human_SINE_Alu_1e-8_8.txt";
my $normalize = $seqFile; 
my $suppressPrint = 0;
my $revcom = 0;
my $range=1; 
my $reverse_editing=0; 
my $hash_ref = analysisSubs::getTripletsAllSeqs($seqFile, $acgt, $siteListFile, $range, $normalize, $reverse_editing, $suppressPrint, $revcom);

exit; 
#calc enrichment scores for best pairs in 1e-8 8 - normalizatoin from within borders of editing sites!
my $freqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-8_8_bestASources/logo_A_Human_SINE_Alu_1e-8_8_freq.txt"; 
my $seqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-8_8_bestASources/seqFasta_A_Human_SINE_Alu_1e-8_8.fa";
my $siteListFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-8_8_bestASources/siteList_A_Human_SINE_Alu_1e-8_8.txt";
my $hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $seqFile, $siteListFile);
# my $hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $seqFile); #no border!!!
analysisSubs::printFreqHash($hash_ref);

exit; 
#borders for triplets normed by BG: 1e-5 5  
my $seqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-8_8_bestASources/seqFasta_A_Human_SINE_Alu_1e-8_8.fa";
my $acgt = "a"; 
my $siteListFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-8_8_bestASources/siteList_A_Human_SINE_Alu_1e-8_8.txt";
my $normalize = $seqFile; 
my $suppressPrint = 0;
my $revcom = 0;
my $range=1; 
my $reverse_editing=0; 
my $hash_ref = analysisSubs::getTripletsAllSeqs($seqFile, $acgt, $siteListFile, $range, $normalize, $reverse_editing, $suppressPrint, $revcom);

exit; 
#borders for triplets normed by BG: 1e-5 5  
my $seqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-5_5/seqFasta_A_Human_SINE_Alu_1e-5_5.fa";
my $acgt = "a"; 
my $siteListFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-5_5/siteList_A_Human_SINE_Alu_1e-5_5.txt"; 
my $normalize = $seqFile; 
my $suppressPrint = 0;
my $revcom = 0;
my $range=1; 
my $reverse_editing=0; 
my $border=1;
my $hash_ref = analysisSubs::getTripletsAllSeqs($seqFile, $acgt, $siteListFile, $range, $normalize, $reverse_editing, $suppressPrint, $revcom, $border);



exit;  
#borders for triplets normed by BG: 1e-8 8 
my $seqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-8_8/seqFasta_A_Human_SINE_Alu_1e-8_8.fa";
# my $seqFile = "tmp.fa";
my $acgt = "a"; 
my $siteListFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-8_8/siteList_A_Human_SINE_Alu_1e-8_8.txt"; 
# my $siteListFile = "siteList.txt"; 
my $normalize = $seqFile; 
my $suppressPrint = 0;
my $revcom = 0;
my $range=1; 
my $reverse_editing=0; 
my $border=1;
# my $hash_ref = analysisSubs::getTriplets($sequence, $acgt, $site_ref, $range, $reverse_editing);
my $hash_ref = analysisSubs::getTripletsAllSeqs($seqFile, $acgt, $siteListFile, $range, $normalize, $reverse_editing, $suppressPrint, $revcom, $border);

exit;  
#check getTerminalsFromSiteList
my $siteListFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-8_8_bestASources/siteList_A_Human_SINE_Alu_1e-8_8.txt";
my $tSites  = analysisSubs::getTerminalsFromSiteList($siteListFile);
foreach my $k(keys %$tSites){
	print $k . "\t".${$tSites->{$k}}[0] . "\t". ${$tSites->{$k}}[1]. "\n"; 
}

exit; 
#calc enrichment scores for best pairs in 1e-8 8 
my $freqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-8_8_bestASources/logo_A_Human_SINE_Alu_1e-8_8_freq.txt"; 
my $seqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-8_8_bestASources/seqFasta_A_Human_SINE_Alu_1e-8_8.fa";
my $hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $seqFile); 
analysisSubs::printFreqHash($hash_ref);

# Result: 
# Position        A       C       G       T
# -2      0.625564162310808       1.2185023068733 1.46193617400606        1.03141115799344
# -1      0.561167227833894       1.91382176192303        0.641895450349141       1.03692672568325
# 0       1       0       0       0
# 1       0.68987341772152        1.25593948378758        1.15154711673699        1.2602855022869
# 2       0.694315803917532       0.939298747808198       1.95401520792437        0.776754890678943

exit; 
 #get best A-sources for Alu 1e-5 5 
my $dir = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-5_5/"; 
my $pairFile = $dir ."graph2_Human_SINE_Alu_1e-5_5.txt"; 
my $fastaS = $dir . "seqFasta_A_Human_SINE_Alu_1e-5_5.fa"; 
my $fastaT = $dir . "seqFasta_G_Human_SINE_Alu_1e-5_5.fa"; 
my $STreversed = 1; 
my $transform = "ga"; 
my $hash_ref = analysisSubs::getBestSources($pairFile, $fastaS, $fastaT, $STreversed, $transform); 

exit; 
 
my $dir = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-8_8/"; 
my $pairFile = $dir ."graph2_Human_SINE_Alu_1e-8_8.txt"; 
my $fastaS = $dir . "seqFasta_A_Human_SINE_Alu_1e-8_8.fa"; 
my $fastaT = $dir . "seqFasta_G_Human_SINE_Alu_1e-8_8.fa"; 
my $STreversed = 1; 
my $transform = "ga"; 
my $hash_ref = analysisSubs::getBestSources($pairFile, $fastaS, $fastaT, $STreversed, $transform); 

exit; 
#check transform to 3-base
my $seqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-8_8/seqFasta_A_Human_SINE_Alu_1e-8_8.fa";
analysisSubs::transformFileTo3base($seqFile, "GA", "transformFileTo3base.output.txt"); 

my $sequence="aaaccctttgggactgactg"; 
my $acgt = "t"; 
my @sites = qw/4 5 6 14 18/; 
my $site_ref = \@sites; 
my $range=1; 
my $reverse_editing=0;

# my $seqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-5_5/seqFasta_A_Human_SINE_Alu_1e-5_5.fa"; 
# my $siteListFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-5_5/siteList_A_Human_SINE_Alu_1e-5_5.txt"; 
my $seqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-8_8_control/seqFasta_A_Human_SINE_Alu_1e-8_8_control.fa";
my $siteListFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-8_8_control/siteList_A_Human_SINE_Alu_1e-8_8_control.txt"; 
my $normalize = $seqFile; 
my $suppressPrint = 0;
my $revcom = 1;
$normalize = 0;
# my $hash_ref = analysisSubs::getTriplets($sequence, $acgt, $site_ref, $range, $reverse_editing);
my $hash_ref = analysisSubs::getTripletsAllSeqs($seqFile, $acgt, $siteListFile, $range, $normalize, $reverse_editing, $suppressPrint, $revcom);


# foreach my $t (sort {$hash_ref->{$a} cmp $hash_ref->{$b}} keys %$hash_ref){
	# print $t ."\t". $hash_ref{$t} ."\n"; 
# }


# analysisSubs::check3("/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-5_5/seqFasta_G_Human_SINE_Alu_1e-5_5.fa");  
# analysisSubs::check3("");  
# analysisSubs::check3("/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-5_5_old-new_1-0/seqFasta_G_Human_SINE_Alu_1e-5_5.fa");  

#(my $hash_ref, my $range )= analysisSubs::readFreqFile("/home/alu/binknis/Data/Human/SINE/results/Tracks/Filtered/tracks_Human_SINE_Alu_1e-5_5_old-new_1-0/logo_A_Human_SINE_Alu_1e-5_5_freq.txt");
#analysisSubs::printFreqHash($hash_ref, "out.temp.txt"); 

#Check getNucFrequencyPerPos function: 
#analysisSubs::check(); 

exit; 
#code for normed motifs in Lizard (for Nika) 
my $freqFile = "/home/alu/binknis/_Data2/LizardNika/LINE/results/Tracks/tracks_Lizard_LINE_L2_1e-8_8/logo_A_Lizard_LINE_L2_1e-8_8_freq.txt"; 
my $seqFile = "/home/alu/binknis/_Data2/LizardNika/LINE/results/Tracks/tracks_Lizard_LINE_L2_1e-8_8/seqFasta_A_Lizard_LINE_L2_1e-8_8.fa"; 
my $hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $seqFile); 
analysisSubs::printFreqHash($hash_ref);

exit; 
#code for finding normalized Alu motifs - for CONTROL
my $freqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-8_8_control/logo_A_Human_SINE_Alu_1e-8_8_control_freq.txt"; 
my $seqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-8_8_control/seqFasta_A_Human_SINE_Alu_1e-8_8_control.fa"; 
my $hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $seqFile); 
analysisSubs::printFreqHash($hash_ref);


exit; 
print "Not trimmed frequencies:\n"; 
my $seqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-8_8/seqFasta_A_Human_SINE_Alu_1e-8_8.fa"; 
my $hash_ref = analysisSubs::getNucFreqFromFasta($seqFile, 2, 0); 
analysisSubs::printFreqHash(\%{$hash_ref->{"a"}});
print "Trimmed frequencies:\n"; 
$seqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-8_8/seqFastaTrimmed_A_Human_SINE_Alu_1e-8_8.fa"; 
$hash_ref = analysisSubs::getNucFreqFromFasta($seqFile, 2, 0); 
analysisSubs::printFreqHash(\%{$hash_ref->{"a"}});

exit; 
#code for finding normalized Alu motif with Trimmed poly-A
my $freqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-8_8/logo_A_Human_SINE_Alu_1e-8_8_freq.txt"; 
my $seqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-8_8/seqFastaTrimmed_A_Human_SINE_Alu_1e-8_8.fa"; 
my $hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $seqFile); 
analysisSubs::printFreqHash($hash_ref);

exit; 
#code for finding normalized Alu motif with Trimmed poly-A
my $freqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-5_5/logo_A_Human_SINE_Alu_1e-5_5_freq.txt"; 
my $seqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-5_5/seqFastaTrimmed_A_Human_SINE_Alu_1e-5_5.fa"; 
my $hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $seqFile); 
analysisSubs::printFreqHash($hash_ref);

#
exit; 
#code for finding normalized Alu motifs
my $freqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-5_5/logo_G_Human_SINE_Alu_1e-5_5_freq.txt"; 
my $seqFile = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-5_5/seqFasta_G_Human_SINE_Alu_1e-5_5.fa"; 
my $hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $seqFile); 
analysisSubs::printFreqHash($hash_ref);


#
exit; 
my $siteList = "/home/alu/binknis/Data/Human/SINE/results/Tracks/tracks_Human_SINE_Alu_1e-5_5/siteList_A_Human_SINE_Alu_1e-5_5.txt"; 
analysisSubs::getTerminalsFromSiteList($siteList); 



 

 