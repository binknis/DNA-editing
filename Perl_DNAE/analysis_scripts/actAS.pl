use strict; 
use lib $ENV{HOME} . "/Perl_DNAE/analysis_scripts"; 
use analysisSubs; 
use lib $ENV{HOME} . "/Perl_DNAE"; 
use getAll; 

my $flags = shift (@ARGV); 
my @serials = split(',', $flags); 
foreach my $s (@serials){
	if ($s==1){
		# print "getNucAndTrioFreq\n"; 
		getNucAndTrioFreq(@ARGV);
	}
	elsif($s==2){
		# print "posInCons\n"; 
		posInCons(@ARGV); 
	}
	elsif($s==3){
		bestSources(@ARGV); 
	}
	elsif($s==4){
		bestSourcesPolyATrimmed(@ARGV);
	}
	elsif($s==5){
		getNucAndTrioFreqRedundant(@ARGV);
	}
	elsif($s==6){
		getPairsByContext(@ARGV); 
	}
	elsif($s==7){
		analysisSubs::getNucsListInSeqForPairs(@ARGV); 
	}
}

sub getNucAndTrioFreq{ 
	(my $dir) = @_; 
	$dir .= "/"; 
	my $GA = "A"; 
	my $ga = lc $GA; 
	my $trim = 0; #0
	my $border = 0; #0
	my $nucleotides = 0; #"cgt"; 
	
	my $seqFile = $dir . "seqFasta_".$GA."_Human_SINE_Alu_1e-0_4.fa";
	if ($trim){
		my $trimmedSeqFile = $seqFile; 
		$trimmedSeqFile =~ s/seqFasta/seqFastaTrimmed/; 
		my $mismatches = 2; #1 
		analysisSubs::trimPolyA($seqFile, $trimmedSeqFile, 0, $mismatches);
		$seqFile = $trimmedSeqFile; 
	}
	
	my $siteListFile = $dir . "siteList_".$GA."_Human_SINE_Alu_1e-0_4.txt"; 
	my $freqFile = $dir . "logo_".$GA.($border ? "_bordered" : "").($nucleotides ? "_".$nucleotides : "").($trim ? "_trimmed" : "")."_Human_SINE_Alu_1e-0_4_freq.txt"; 
	my $range = 2; 
	my $cluster_file = 0; #$dir ."clusters_Human_SINE_Alu_1e-0_4.tab|1"; #***0
	if ($cluster_file){
		analysisSubs::getNucFrequencyPerPosAllSeqs($seqFile, $cluster_file, $range, 1, 0, $freqFile); #(my $seqFile, my $siteListFile, my $range, my $normalize, my $reverse_editing, my $outfile)
	}
	else{
		analysisSubs::getNucFrequencyPerPosAllSeqs($seqFile, $siteListFile, $range, 1, 0, $freqFile); #(my $seqFile, my $siteListFile, my $range, my $normalize, my $reverse_editing, my $outfile)
	}
	my $hash_ref; 
	if ($border){
		$hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $seqFile, $siteListFile, $nucleotides); #norm within borders!!!
	}
	else{
		$hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $seqFile, 0, $nucleotides);
	}
	my $freqNormedFile = $dir . "freq_".$GA.($border ? "_bordered" : "").($nucleotides ? "_".$nucleotides : "").($trim ? "_trimmed" : "").($cluster_file ? "_redundant" : "")."_Human_SINE_Alu_1e-0_4.txt"; 
	analysisSubs::printFreqHash($hash_ref, $freqNormedFile);

	my $acgt = $ga; my $normalize = $seqFile; my $suppressPrint = 1; my $revcom = 0; my $range=1; my $reverse_editing=0;  
	# my $acgt = $ga; my $normalize = $seqFile; my $suppressPrint = 0; my $revcom = 0; my $range=1; my $reverse_editing=0; my $border=1; #norm within borders!!!
	my $triplets_ref = analysisSubs::getTripletsAllSeqs($seqFile, $acgt, $siteListFile, $range, $normalize, $reverse_editing, $suppressPrint, $revcom, $border);
}

sub getNucAndTrioFreqRedundant{
	(my $dir) = @_; 
	$dir .= "/"; 
	my $seqFile = $dir . "seqFasta_A_Human_SINE_Alu_1e-0_4.fa";

	my $siteListFile = $dir . "siteList_A_Human_SINE_Alu_1e-0_4.txt"; 
	my $clustersFile = $dir . "clusters_Human_SINE_Alu_1e-0_4.tab"; 
	my $siteListOrClusterFile = $clustersFile."|1";   #1 - for A sites, 0 - for G sites
	my $freqFileRedundant = $dir . "freqRedundant_A_Human_SINE_Alu_1e-0_4.txt"; 
	my $range = 2; 
	analysisSubs::getNucFrequencyPerPosAllSeqs($seqFile, $siteListOrClusterFile, $range, 1, 0, $freqFileRedundant); #(my $seqFile, my $siteListFile, my $range, my $normalize, my $reverse_editing, my $outfile)
	
	my $hash_ref = analysisSubs::normalizeFreqByFasta($freqFileRedundant, $seqFile); 
	# my $hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $seqFile, $siteListFile); #norm within borders!!!
	analysisSubs::printFreqHash($hash_ref);

#	my $acgt = "a"; my $normalize = $seqFile; my $suppressPrint = 1; my $revcom = 0; my $range=1; my $reverse_editing=0; my $border=0; 
#	# my $acgt = "a"; my $normalize = $seqFile; my $suppressPrint = 0; my $revcom = 0; my $range=1; my $reverse_editing=0; my $border=1; #norm within borders!!!
#	analysisSubs::getTripletsAllSeqs($seqFile, $acgt, $siteListOrClusterFile, $range, $normalize, $reverse_editing, $suppressPrint, $revcom, $border);
}

sub posInCons{
	(my $dir) = @_; 
	$dir .= "/"; 
	my $siteListFile = $dir ."/". "siteList_A_Human_SINE_Alu_1e-0_4.txt";
	analysisSubs::getEditedPositionsInCons($siteListFile);
	my $nucListFile = $dir ."/". "nucList_A_Human_SINE_Alu_1e-0_4.txt";
	analysisSubs::nucListToFreq($nucListFile); 
}

sub bestSourcesPolyATrimmed {
	(my $dir) = @_; 
	$dir .= "/"; 
	#trim source-file
	my $seqFileS = $dir . "seqFasta_A_Human_SINE_Alu_1e-0_4.fa";
	my $trimmedSeqFileS = $dir . "seqFastaTrimmed_A_Human_SINE_Alu_1e-0_4.fa";
	my $mms = 2; 
	analysisSubs::trimPolyA($seqFileS, $trimmedSeqFileS, 0, $mms); 
	#trim target-file
	my $seqFileT = $dir . "seqFasta_G_Human_SINE_Alu_1e-0_4.fa";
	my $trimmedSeqFileT = $dir . "seqFastaTrimmed_G_Human_SINE_Alu_1e-0_4.fa";
	analysisSubs::trimPolyA($seqFileT, $trimmedSeqFileT, 0, $mms); 
	#get best sources from trimmed file 
	my $pairFile = $dir ."graph2_Human_SINE_Alu_1e-0_4.txt";
	my $fastaS = $trimmedSeqFileS;
	my $fastaT = $trimmedSeqFileT;
	my $STreversed = 1; #for A>I 
	# my $STreversed = 0; #for G>A 
	my $transform = "ga"; 
	my $trimmed = 1; 
	my $hash_ref = analysisSubs::getBestSources($pairFile, $fastaS, $fastaT, $STreversed, $transform, $trimmed);
}

sub bestSources {
	(my $dir) = @_; 
	$dir .= "/"; 
	my $suffix = "Lizard_LINE_L1_1e-8_8"; 
	#my $suffix = "Human_SINE_Alu_1e-0_4"; 
	my $pairFile = $dir ."graph2_".$suffix.".txt";
	#my $fastaS = $dir . "seqFasta_A_".$suffix.".fa"; #A>I
	#my $fastaT = $dir . "seqFasta_G_".$suffix.".fa"; #A>I
	my $fastaS = $dir . "seqFasta_G_".$suffix.".fa"; #G>A
	my $fastaT = $dir . "seqFasta_A_".$suffix.".fa"; #G>A
	#my $STreversed = 1; #A>I
	my $STreversed = 0; #G>A
	#my $transform = "ga"; 
	my $transform = 0; 
	my $trimmed = 0; 
	my $hash_ref = analysisSubs::getBestSources($pairFile, $fastaS, $fastaT, $STreversed, $transform, $trimmed);
}


sub getPairsByContext {
	(my $dir) = @_; 
	$dir .= "/"; 	
	# my $clusterTabFile = $dir ."clusters_Human_SINE_Alu_1e-0_4".($dir =~ /_control/ ? "_control" : "").".tab"; 
	# my $SorT = 1; #0- Gseq, 1- Aseq
	# my $position = -1; 
	my $clusterTabFile = $dir ."clusters_Lizard_LINE_L1_1e-8_8".($dir =~ /_control/ ? "_control" : "").".tab"; 
	my $SorT = 1; #0- Gseq, 1- Aseq
	my $position = -1; 
	analysisSubs::getPairsByContext($clusterTabFile, $SorT, $position);
	
}




# if ($p = getNucFrequencyPerPos)
# if ($p = getNucFrequencyPerPosAllSeqs
# if ($p = getNucFreqFromFasta
# if ($p = check3
# if ($p = normalizeFreqByFasta
# if ($p = readFreqFile
# if ($p = printFreqHash
# if ($p = trimPolyA
# if ($p = countSeqSitesInTrimmed
# if ($p = getTriplets
# if ($p = getTripletsAllSeqs
# if ($p = revcomHashKeys
# if ($p = revcom
# if ($p = getTripletsBG
# if ($p = revEditing
# if ($p = check
# if ($p = alignment_to_mm_arr
# if ($p = check2
# if ($p = motif_per_alignment
# if ($p = sitesFromSiteList
# if ($p = nucsFromNucsList
# if ($p = getTerminalsFromSiteList
# if ($p = transformFileTo3base
# if ($p = transformSeqTo3base
# if ($p = getBestSources
# if ($p = getPairs
# if ($p = getEditedNucsInCons
# if ($p = getEditedPositionsInCons
# if ($p = splitFastaBySubfam
# if ($p = getConsPosHist
# if ($p = idFromSubfam
# if ($p = writeSubfamHist
# if ($p = nucListToFreq
# if ($p = markEditingSitesInBLAST
# if ($p = redundantSiteList
# if ($p = getPairSites
# if ($p = uniq
# if ($p = createConsensus
# if ($p = getOldVsNew



# my @ar = qw/1 2 3 4/; 
# print @ar , "\n"; 
# my $a = shift @ar; 

# print @ar, "\n"; 
# print $a ."\n"; 




