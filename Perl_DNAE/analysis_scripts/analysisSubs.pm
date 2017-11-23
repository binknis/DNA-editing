use strict; 
package analysisSubs; 
use Math::Round qw(nearest);
use Bio::SearchIO;
use Bio::SeqIO; 

use List::Util qw(max min sum);
use File::Path qw(mkpath);
use File::Copy;
# use Data::Dumper;

# use lib $ENV{HOME} . "/Perl_DNAE"; #*** make sure it works after replacing this with following lines
use FindBin;
use lib "$FindBin::Bin/..";  # use the parent directory of analysis_scripts directory from where this is run
use lib "$FindBin::Bin"; #because this file is in analysis_scripts, this includes that directory
use getAll; 

###Specifications: 
# splitFastaBySubfam(), getConsPosHist() - need deflines to have subfam at end of defline, with a preceding '=', so that matches: /=([^=]+)$/ to extract subfam


###*** To do: 
#1. check no_CpG filter in getNucFrequencyPerPos
#2. All scripts that read tabClustFile need to make sure that if filter is activated that they only read clusters that are in filter


### Subroutine getNucFrequencyPerPos ###
#Function: gets the nucleotide frequency in positions in relation to a list of sites in a given sequence
#Input: $sequence = a sequence (lc/uc)
#		$site_ref = an array of sites in the sequence to screen positions related to them (1-base)
#		$range = the position range to check (-range .. range) except for 0
#		$normalize (flag)= 	0 - no normalization (returns count); 
#							1 - return fraction of each nuc in each pos; 
#		$reverse_editing = option to reverse_editing (change edited positions back to pre-editing form; e.g. "GA" will reverse A to G); 
#		$no_CpG = Don't count sites that have specific nucs in their proximity(CpG is just an example and typical usage)
#				examples: For no C at -1 arg: "-1,c"; For no T at +2 arg: "1,t"
#Returns: a hash with the frequencies in hash{$pos}{$nuc} = $freq
#Note: 1. if normalizes then rounded to 3 decimal digits
#
#General note: Tabular input: [0]mmType, [1]assembly, [2]class, [3]family, [4]subfam, [5]coordsS, [6]coordsT, [7]num_edited_sites, [8]num_all_mms, [9]total_prob, [10]whereS, [11]whereT, [12]num_clusts, [13]mmSerials, [14]clusters_span_woGaps, [15]clusters_span, [16]alignment_len, [17]mmCount_str
#*** check if no_CpG filter works; the rest of the code isn't affected
sub getNucFrequencyPerPos{
	(my $sequence, my $site_ref, my $range, my $normalize, my $reverse_editing, my $no_CpG) = @_;
	my $noPos, my $noNuc; 
	if ($no_CpG){
		($noPos, $noNuc) = split (',',$no_CpG); 
		$noNuc = lc $noNuc; 
	}
	
	#convert sequence to array and modify it
	my $seq; 
	if ($reverse_editing){ #change edited nucs back to pre-editing nuc
		$sequence = revEditing($sequence, $site_ref, $reverse_editing); 
	}
	else{ #don't reverse editing
		$sequence = lc $sequence;
	}
	my @seq = split('',$sequence);
	
	#init variables
	my @nucs = ('a', 'c', 'g', 't'); 
	my @positions = (-$range .. $range);
	my %nucFreq = (); 
	my %posCount = (); #counter for number of positions counted for each position (isn't necessarily equal to the number of sites because some sites are close to start or end of sequence and won't have some of the positions)
	foreach my $pos (@positions){
		foreach my $nuc (@nucs){
			$nucFreq{$pos}{$nuc}=0;
		}
		$posCount{$pos}=0;
	}
	
	#get frequencies
	foreach my $site(@$site_ref){
		if ($no_CpG){ #check no-CpG filter to filter out specific sites (CpG is just an example and typical usage)
			my $noInd = $site + $noPos - 1; 
			if($noInd <= $#seq && $noInd > 0){ #in bounds of sequence
				next if ($seq[$noInd]) eq $noNuc; #nuc to skip detected
			}
		}
		foreach my $pos (@positions){
			my $ind = $site + $pos - 1; 
			if ($ind <= $#seq && $ind > 0){	#in bounds of sequence
				if ($seq[$ind] =~ /[acgt]/){ #not n or some other character
					$nucFreq{$pos}{$seq[$ind]}++;
					$posCount{$pos}++; 
				}
			}
		}
	}

	#normalize
	if ($normalize > 0){ 
		foreach my $pos (@positions){
			foreach my $nuc (@nucs){
				$nucFreq{$pos}{$nuc} /= $posCount{$pos} if $posCount{$pos};
				$nucFreq{$pos}{$nuc} = nearest(.001, $nucFreq{$pos}{$nuc}); #round to max of 3 digits after the decimal number
			}
		}
	}
	return (\%nucFreq);
}

#Function: sum nuc frequencies for all seqs in fasta file in relation to specific (editing) positions
#Input: 1. seqFile = fasta file for which to calc nuc freq
#		2. siteListFile = 	opt1: siteList file with editing sites for every sequence in the file
#		   					opt2: tabular cluster file and the 1/0 flag for target/source 
#		3. range = How many bases on each side of 'center' nuc to calculate
#		4. normalize = flag: 0- return frequencies, 1- return fractions of each nuc in each pos
#		5. reverse_editing = if to reverse editing sites (e.g. "ga" will reverse 'a' back to 'g')
#		5. outfile = output file for frequencies (if not specified - prints to stdout; if "none" - doesn't print to file)
#		6. no_CpG = Don't count sites that have specific nucs in their proximity(CpG is just an example and typical usage)
#		7. freqPerSeqFile = File to print frequencies for every sequence to (will print only if specified)
#Note: Fasta file and siteList file deflines must be identical!
sub getNucFrequencyPerPosAllSeqs{
	(my $seqFile, my $siteListFile, my $range, my $normalize, my $reverse_editing, my $outfile, my $no_CpG, my $freqPerSeqFile) = @_;
	$range = 2 unless $range; #default - 2 nucs on each side
	my %freqAllSeqs = ();
	my %freqPerSeqAll = (); 
	my %freqCount = ();
	my @poses = (-$range .. $range);
	#Get editing sites from file
	my $sites_r; 
	my $redundant = ($siteListFile =~ /clusters_/ ? 1 : 0); 
	if ($redundant){ #get redundant sites from tabular cluster file
		(my $cluster_file, my $SorT) = split ('\|', $siteListFile); 
		$sites_r = redundantSiteList($cluster_file, $SorT); 
	}
	else{ #get non-redundant sites from siteList file
		$sites_r = sitesFromSiteList($siteListFile); 
	}
	
	# ***temp print: 
		# print "cat siteListFile: $siteListFile\n"; 
		# system("cat $siteListFile"); 
		# print "print siteListFile: $siteListFile\n"; 
		# foreach my $s (keys %$sites_r){
			# print $s, " ", join (',', @{$sites_r->{$s}}) , "\n"; 
		# }
	# *** end temp 
	
	#Get nucFreq around editing sites of each seq in file
	my $inseq = Bio::SeqIO->new( -file => $seqFile, -format => 'fasta' );
	while ( my $seq = $inseq->next_seq ){
		#get freq per seq
		my $s = lc $seq->seq(); 
		my $id; 
		($id) = $seq->id =~ /([^=]+:\d+-\d+[+-])/; #IDs in hash are coords, so get them from seq defline
		# print $id . "doesn't exist\n" unless exists $sites_r->{$id};
		next unless exists $sites_r->{$id}; #skip sequences that aren't in editing list.
		# print $id . "exists\n"; 
		my $freqPerSeq = getNucFrequencyPerPos($s, $sites_r->{$id}, $range, 0, $reverse_editing); 
		$freqPerSeqAll{$id} = $freqPerSeq; #save count per each sequence for printing to per-seq file
		# printFreqHash($freqPerSeq); #***temp
		#sum frequencies 
		foreach my $p (keys %$freqPerSeq){ #each position
			foreach my $n (keys %{$freqPerSeq->{$p}}){ #each nucleotide
				$freqAllSeqs{$p}{$n} += $freqPerSeq->{$p}{$n}; 
				$freqCount{$p} += $freqPerSeq->{$p}{$n}; #sum for normalization (for fraction; not used for enrichment)
			}
		}
	}
	## Print freq for each sequence
	if ($freqPerSeqFile){
		printFreqPerSeq($freqPerSeqFile, \%freqPerSeqAll);
	}
	
	## normalize ##
	if ($normalize){ #create fraction
		foreach my $p (keys %freqAllSeqs){ #each position
			foreach my $n (keys %{$freqAllSeqs{$p}}){ #each nucleotide
				$freqAllSeqs{$p}{$n} /= $freqCount{$p} if $freqCount{$p};
				$freqAllSeqs{$p}{$n} = nearest(.001, $freqAllSeqs{$p}{$n}); #round to max of 3 digits after the decimal number
			}
		}
	}
	
	#print frequencies ($outfile not specified) or write to file ($outfile specified) and return as hash
	printFreqHash(\%freqAllSeqs, $outfile) unless $outfile eq "none";
	# printFreqHash(\%freqAllSeqs); #***temp instead of above
	return \%freqAllSeqs; 
}


#Function: get the frequency of each nucleotide in relation to a specific letter
#		Primary usage is for calculating background frequencies for motif inference of editing sites
#Input: 1. Fasta file (from which to calculate)
#		2. Range: amount of sites before and after each 
#		3. retfmt: return value format - 0 or "" (default): frequency; 
#										 1 = fraction of each nuc in the relative position 
#		4. siteList: the list of edited sites to border the regions from which to calculate frequency (from 1st edited to last edited + flanks of length 'range')
#RetVal: a hash with 3 levels: 	1. {a,c,t,g} - the nucleotide in center position
#								2. {pos} - number from -range to +range
#								3. {a,c,t,g} - each nucleotide (the nucs in each position)
#								4. Value =  frequency, fraction or percent of each nuc in the position (rounded to 3-digits after decimal)
sub getNucFreqFromFasta {
	(my $seqFile, my $range, my $retfmt, my $siteList) = @_;
	if (not -e $seqFile){ #file doesn't exist
		print "$seqFile doesn't exist\n"; 
		return 0; 
	}
	else{
		#print "else in getNucFreqFromFasta\n"; 
	}	
	### get nucleotide frequencies in each position ###
	my @nucs = ('a', 'c', 'g', 't'); 
	my @positions = (-$range .. $range);
	my %nucFreq = ();
	my %sum = (); #for converting frequecies to fractions
	foreach my $center (@nucs){ #center is nuc of reference which we seek its surrounding motif
		foreach my $pos (@positions){
			$sum{$center}{$pos}=0; 
			foreach my $nuc (@nucs){
				$nucFreq{$center}{$pos}{$nuc}=0;
			}
		}
	}
	
	### get border per sequence, based on editing sites ### 
	my $borders_r = ""; 
	if ($siteList){
		$borders_r = analysisSubs::getTerminalsFromSiteList($siteList);
	}
	
	### fetch nucleotide frequencies around each nucleotide ###
	my $inseq = Bio::SeqIO->new( -file => $seqFile, -format => 'fasta' );
	while ( my $seq = $inseq->next_seq ){
		my $seqToSplit = lc $seq->seq();
		if ($borders_r){ #get substr within first and last editing site -- in order to calc background from seq within borders of first and last editing sites (with immediate flanks of length $range)
			(my $coords) = $seq->id() =~ /([^=|]+:\d+-\d+[+-])/;
			my $flank = $range; #seq used should include the flanks immediately next to first and last editing site. 
			my $start = $borders_r->{$coords}[0] - $flank; #include left flank of 1st editing site. 
			$start = 0 if $start < 0;
			my $end = $borders_r->{$coords}[1]; 
			my $substrLen = $end - $start + 1 + $flank; #substr length should include end site and right flank of last editing site. 
			$seqToSplit = substr($seqToSplit, $start, $substrLen); #(note: $start + $substrLen may be out-of-bounds - in that case substr returns only until the last char)
		}
		my @s = split("", $seqToSplit);
		for (my $i=0; $i<=$#s; $i++){ #i is the nuc of reference in center of motif
			next if $s[$i] !~ /[actg]/; #skip n/N etc.
			for (my $j=max(0,$i-$range); $j<=min($i+$range, $#s); $j++){
				next if $s[$j] !~ /[actg]/; 
				$nucFreq{$s[$i]}{$j-$i}{$s[$j]}++; #{center-nuc}{position relative to center}{nuc found}
			}
		}
	}
	
	#Alter output format(optional)
	if ($retfmt){
		#sum nuc freqs in each pos
		foreach my $cNuc (@nucs){
			foreach my $pos (@positions){
				foreach my $nuc (@nucs){
					$sum{$cNuc}{$pos} += $nucFreq{$cNuc}{$pos}{$nuc}; 
				}
			}
		}
		#convert from freqs to fractions
		foreach my $cNuc (@nucs){ 
			foreach my $pos (@positions){
				foreach my $nuc (@nucs){
					next if $sum{$cNuc}{$pos}==0; #avoid division by 0
					$nucFreq{$cNuc}{$pos}{$nuc} /= $sum{$cNuc}{$pos};
					$nucFreq{$cNuc}{$pos}{$nuc} = nearest(.001, $nucFreq{$cNuc}{$pos}{$nuc}); #round to max of 3 digits after the decimal number	
				}
			}
		}
		
	}
	return \%nucFreq; 
}

#Function: normalizes one nuc-frequency for positions relative to a certain nuc by another nuc frequency
#input: 1. file with motif frequencies from editing data
#		2. fasta file used to extract the background frequencies for normalization
#		3. siteList = use borders 
#		4. nucleotides = which nucleotides to use for output (e.g. "cgt" for C, G and T only. good to filter biases caused by polyA)
#			Note: if "nucleotides" doesn't include the center-Nuc (the edited nuc; see below) then "a" will be considered the cNuc 
#retval: hash-ref with normalized freqs; hash format: HASH{pos-relative-to-center-of-motif-nuc}{nuc}=freq; e.g. HASH{-1}{'c'}=0.432
sub normalizeFreqByFasta{
	(my $freqFile, my $seqFile, my $siteListFile, my $nucleotides) = @_;
	#Extract frequencies from freq-file
	(my $freq_r, my $range) = readFreqFile($freqFile); 
	die "File not found: $freqFile\n" unless $freq_r; 
	
	#Retain only nucs of interest in nuc-array and create poses array
	my @nucs = ('a', 'c', 'g', 't');
	$nucleotides = lc $nucleotides; 
	@nucs = split(//, $nucleotides) if $nucleotides; #to get frequencies for only some nucleotides (see above)
	my @poses = (-$range .. $range);
	
	##get background frequencies for normalization
	my $bgFreq_r = getNucFreqFromFasta($seqFile, $range, 1, $siteListFile); 
	###*** Temp - print to see nucFreq
	# foreach my $center ("a", "g"){ 
		# print $seqFile ."\n". "center nuc: $center\n"; 
		# printFreqHash($bgFreq_r->{$center});
	# }
	###*** END Temp
	#print "bgFreq_r:\n"; #***
	#printFreqHash($bgFreq_r->{"a"}); #***
	#extract the nuc in center of motif from the freq-file's data (to avoid sending it as parameter; e.g. if the line is "0 0 1 0" then the center nuc is G)
	my $cNuc; 
	for (my $i=0; $i<=$#nucs; $i++){
		# print $nucs[$i] . "\n"; 
		$cNuc = $nucs[$i] if $freq_r->{"0"}{$nucs[$i]}==1;  #"0" refers to position in center of motif (see "readFreqFile")
	}
	$cNuc = "a" if ($nucleotides and not $cNuc); #ONLY CORRECT WHEN "a" IS REMOVED IN "NUCLEOTIDES" LIST
	die "No center-nuc found for $freqFile - exiting normalizeFreqByFasta\n" unless defined $cNuc;

	#Retain only nucleotides of interest in freq_r 
	if ($nucleotides){ 
		#delete unwanted
		foreach my $nuc ('a', 'c' , 'g', 't'){
			unless ($nucleotides =~ /$nuc/){ #nucleotide wasn't specified
				foreach my $pos (@poses){ #del from all positions in hashRef
					delete $freq_r->{$pos}{$nuc}; 
				}
			}
		}
		
		#rescale remaining fracs to sum-1 per pos
		my %posSum = ();
		my %posSum_bg = ();
		foreach my $pos (@poses){ 
			$posSum{$pos}=0;
			$posSum_bg{$pos}=0;
			foreach my $nuc (@nucs){ #sum
				$posSum{$pos} += $freq_r->{$pos}{$nuc};
				$posSum_bg{$pos} += $bgFreq_r->{$cNuc}{$pos}{$nuc}; 
			}
			foreach my $nuc (@nucs){ #rescale to 1-sum
				$freq_r->{$pos}{$nuc} /= $posSum{$pos} if $posSum{$pos};
				$bgFreq_r->{$cNuc}{$pos}{$nuc} /= $posSum_bg{$pos} if $posSum_bg{$pos};
			}
		}
	
	}
	
	## Normalize (generate enrichment score)
	foreach my $pos (@poses){ 
		foreach my $nuc (@nucs){
			if ($bgFreq_r->{$cNuc}{$pos}{$nuc}){ #divide by background
				$freq_r->{$pos}{$nuc} /= $bgFreq_r->{$cNuc}{$pos}{$nuc};
				$freq_r->{$pos}{$nuc} = nearest(0.001, $freq_r->{$pos}{$nuc}); 
			}
			elsif ($freq_r->{$pos}{$nuc}>0 and $bgFreq_r->{$cNuc}{$pos}{$nuc}==0) {#avoid dividing by zero
				$freq_r->{$pos}{$nuc} = "inf"; #not middle line, where we expect 0
			}
			else{ #the freq is 0 anyway, but for completeness I added the command
				$freq_r->{$pos}{$nuc}=0; 
			}
		}
	}
	
	return $freq_r; 
}

#Function: normalizeFreqfileByFreqfile - get background frequencies (nuc-distributions) for all sites that weren't edited in the edited elements
#Method: Get frequency (count) of all sites in background, deduct frequencies (count) derived from editing sites only (then optionally convert to nuc-distribution)
#*** CHECK
sub getBackgroundNoEditedFreq {
	(my $GA, my $seqFile, my $siteListFile, my $range, my $seqFile2, my $normalize, my $outfile, my $reverse_editing) = @_;
	my $retfmtInSubs = 0; 
	$seqFile2 = $seqFile unless $seqFile2; #default: seqFile1 is like seqFile2
	my %finalFreq = ();
	
	#Get editing sites frequencies 
	my $editedFreq = getNucFrequencyPerPosAllSeqs($seqFile, $siteListFile, $range, $retfmtInSubs, $reverse_editing, "none");  
	# print "editedFreq:\n"; #***TEMP
	# printFreqHash($editedFreq); #***TEMP
	#Get background frequencies - all sites in fasta file
	my $bgFreqAll = getNucFreqFromFasta($seqFile2, $range, $retfmtInSubs); 
	my $bgFreq = $bgFreqAll->{lc $GA}; 
	# print "bgFreq:\n"; #***TEMP
	# printFreqHash($bgFreq); #***TEMP
	
	#subtract editing site frequencies 
	foreach my $p (keys %$bgFreq){ #each position
		foreach my $n (keys %{$bgFreq->{$p}}){ #each nucleotide
			$finalFreq{$p}{$n} = $bgFreq->{$p}{$n} - $editedFreq->{$p}{$n};
		}
	}
	
	## normalize ##
	if ($normalize){ #create fraction
		my %freqCount = (); 
		#sum frequencies and divide by sum for fractions
		foreach my $p (keys %finalFreq){ #each position
			foreach my $n (keys %{$finalFreq{$p}}){ #each nucleotide
				$freqCount{$p} += $finalFreq{$p}{$n}; 
			}
			#divide by sum per position
			foreach my $n (keys %{$finalFreq{$p}}){ #each nucleotide
				$finalFreq{$p}{$n} /= $freqCount{$p} if $freqCount{$p};
				$finalFreq{$p}{$n} = nearest(.001, $finalFreq{$p}{$n}); #round to max of 3 digits after the decimal number
			}
		}
	}
	
	#print frequencies if $outfile was specified($outfile not specified) and return results as hash
	printFreqHash(\%finalFreq, $outfile) if $outfile;
	# printFreqHash(\%finalFreq); #***TEMP
	return \%finalFreq; 
}



#Function(helper): Reads a frequence file and returns the data as a hash (format is as in logo_A_..._freq.txt files)
#Input: the freq file
#RetVal: 1. hash with the freq file data; 
#		 2. range of motif checked (+/- how many bases)
sub readFreqFile{
	(my $freqFile) = @_; 
	my %freqs = (); 
	my @nucs = ('a', 'c', 'g', 't'); 
	open (my $freq_fh, $freqFile) or return 0; 
	while(my $line = <$freq_fh>){
		chomp $line;
		if ($line =~ /^(-?\d+)\t/){ #skip 1st line
			my @f = split(/\t/, $line); #first field is the position in relation to edited site, next 4 are a,c,g,t freqs
			my $pos = shift @f; 
			for(my $i=0; $i<=3; $i++){
				$freqs{$pos}{$nucs[$i]} = $f[$i]; 
			}
		}
	}
	close($freq_fh); 
	my $range = max(keys %freqs); 
	return (\%freqs, $range); 
}

#Function(helper): print a frequence hash to output or file
sub printFreqHash{
	(my $freq_ref, my $outFile) = @_; 
	my $range = max(keys %$freq_ref); 
	my @poses = (-$range .. $range);
	my $fh;
	if ($outFile) {
	   open($fh, '>', $outFile) or die;
	} else {
	   $fh = \*STDOUT;
	}
	print  $fh "Position";
	foreach my $n (sort keys %{$freq_ref->{"0"}}){
		print $fh "\t". uc $n;
	}
	print $fh "\n"; 
	foreach my $pos (@poses){
		print $fh $pos ."\t"; 
		foreach my $nuc ('a','c','g','t'){
			if (exists $freq_ref->{$pos}{$nuc}){ #Should happen only when normalizing the freq hash with only some nucleotides (see $nucleotides in normalizeFreqByFasta)
				print $fh $freq_ref->{$pos}{$nuc} . ($nuc eq 't' ? "\n" : "\t"); 
			}
		}
	}
	close($fh) if $outFile; 
}

#Function: printFreqPerSeq - prints frequency for every sequence to file (helper for getNucFrequencyPerPosAllSeqs)
sub printFreqPerSeq{
	(my $fpsFile, my $fpsAll) = @_; 
	my @nucs = qw/a c g t/; 
	my %fps; 
	#open outfile and print header
	open (my $fps_fh, ">". $fpsFile) || die "open $fpsFile\n"; 
	my $header = "id" . "\t". "pos" ."\t". uc(join("\t", @nucs)); 
	print $fps_fh $header ."\n";
	
	foreach my $id (keys %$fpsAll){ #each seqID
		%fps = %{$fpsAll->{$id}};
		#print frequencies per string		
		foreach my $p (sort {$a <=> $b} keys %fps){ #each position
			my $outStr = $id . "\t" . $p; 
			foreach my $n (@nucs){ #each nucleotide
				$outStr .= "\t" . $fps{$p}{$n}; 
			}
			print $fps_fh $outStr ."\n"; 
		}	
	}
	close($fps_fh); 
}

#Function: trimPolyA - trim polyA from sequences and return original deflines (by default, optionally doesn't retain)
sub trimPolyA {
	(my $seqFile, my $outFile, my $dontRetainDeflines, my $mismatches) = @_; 
	$mismatches = 1 unless $mismatches; #amount of consecutive mismatches in polyA allowed and still continues removing (default = 1, like trimest's default).
	system ("/private/common/Software/EMBOSS/EMBOSS-6.4.0/emboss/bin/trimest -sequence $seqFile -outseq $outFile -mismatches $mismatches"); 
	return if $dontRetainDeflines; 
	#get original deflines
	my $counter = 0; 
	my @ids = (); 
	my $inseq = Bio::SeqIO->new( -file => $seqFile, -format => 'fasta' );
	while ( my $seq = $inseq->next_seq ){
		$ids[$counter] = $seq->id; 
		$counter++;
	}
	#create temp output name (concatenation of input seq-file and trimmed file, in input seq-file's dir)
	my $tempSuffix = $outFile; 
	$tempSuffix =~ s/\S+\///; 
	my $tempFile = "temp.". $inseq . "." . $tempSuffix;
	
	#replace trimmed-file deflines (to leave same as original - needed for other functions)
	my $trimmed = Bio::SeqIO->new( -file => $outFile, -format => 'fasta' );
	my $outseqTemp = Bio::SeqIO->new( -file => ">$tempFile", -format => 'fasta' ); #temp outfile with 
	$counter = 0; 
	while ( my $seq = $trimmed->next_seq ){
		$seq->display_id($ids[$counter]); #change ID
		$outseqTemp->write_seq($seq); #write to temp file
		$counter++; 
	}
	#move temp file back to original output and delete temp file
	# unlink($outfile); 
	File::Copy::move($tempFile, $outFile) || die "Didn't move $tempFile back to $outFile\n"; 
	unlink($tempFile);
}

#Function: countSeqSitesInTrimmed - counts amount of sequences that are present in a modified sequence file
#Why it was written: Because motif was biased towards A in -2..+2, I calculated motif after trimming polyA. 
#	Editing sites in polyA were discarded so I wanted to see which were retained after the trimming. 
#Output: number of sequences with any sites still present and number of 
sub countSeqSitesInTrimmed{
	(my $seqFile, my $siteListFile) = @_; 
	my $seqCount = 0; 
	my $siteCount = 0; 
	my $sites_ref = sitesFromSiteList($siteListFile); 
	my $trimmed = Bio::SeqIO->new( -file => $seqFile, -format => 'fasta' );
	while ( my $seq = $trimmed->next_seq ){
		my $sites=0; 
		(my $coords) = $seq->display_id =~ /([^=|]+:\d+-\d+[+-])/; 
		foreach my $site(@{$sites_ref->{$coords}}){ 
			if ($site <= $seq->length){ #(sites are 1-based)
				$sites++; 
			}
		}
		if ($sites){
			$seqCount++; 
			$siteCount += $sites; 
		}
	}
	print "seqs\tsites\n$seqCount\t$siteCount\n"; 
}


#Function: gets the triplet-motif of nucs around the editing sites in one sequence
#Input: 1. sequence = nuc string
#		2. $acgt = nuc for center of triplet
#		3. site_ref = array with editing sites (get motifs around these sites)(1-base)
#		4. range = 1 for triplets, 2 for quadruplets, etc... 
#		5. reverse_editing = if to reverse editing sites
#Note: 1. triplets is just the default but could be used for any uneven number >=3
#		2. Doesn't normalize because caller function needs to sum (and normalize)
#		3. triplets with any character other than acgt will be skipped
#		4. Possible bug: if reverese_editing then middle nuc in motif may have been reversed.
sub getTriplets{
	(my $sequence, my $acgt, my $site_ref, my $range, my $reverse_editing) = @_;
	$range = 1 unless $range; #default motif size = triplets
	#reverse editing in seq and convert to array
	if ($reverse_editing){ #change edited nucs back to pre-editing nuc (also converts to lc)
		$sequence = revEditing($sequence, $site_ref, $reverse_editing); 
	}
	else{ #don't reverse editing
		$sequence = lc $sequence;
	}
	my @seq = split('',$sequence);
	$_-- for @$site_ref; #convert sites to 0-base for ease of use
	my @tmp = @$site_ref; 
	# print "@tmp\n"; #***
	# get frequencies
	my %triFreq = (); 	
	my $tri; 
	foreach my $site (@$site_ref){
		if ($site-$range >= 0 && $site+$range <= $#seq){ #in bounds of sequence
			$tri = join ('', @seq[$site-$range .. $site+$range]);
			next unless ($tri =~ /^[acgt]+$/); #no n or other characters in triplet
			$triFreq{$tri}++;
		}
	}
	# foreach my $t (sort {$triFreq{$b} <=> $triFreq{$a}} keys %triFreq){
			# print $t ."\t". $triFreq{$t} ."\n"; 
	# }
	return (\%triFreq);
}



#Function: sum triplet frequencies for all seqs in fasta file
#Input: 1. seqFile = fasta file for which to calc triplets
#		2. $acgt = nuc for center of triplet
#		3. siteListFile = 	opt1: siteList file with editing sites for every sequence in the file
#		   					opt2: tabular cluster file and the 1/0 flag for target/source 
#		4. range = 1 for triplets, 2 for quadruplets, etc... 
#		5. normalize = flag: 0- return frequncies, 1- return fractions of all triplets, 2 - fileName - enrichment over triplets in the file (will return enrichment score!)
#		6. supressPrint = suppress printing to stdout
# 		7. $revcom = revcom triplets before returning and/or printing (helpful for comparing motif of sense and anti-sense editing)
#		8. border = flag: 0- calc BG from whole sequences. 1- calc BG for normalization using substr in borders of 1st and last editing site only.
#Note: Fasta file and siteList file deflines must be identical!
sub getTripletsAllSeqs{
	(my $seqFile, my $acgt, my $siteListFile, my $range, my $normalize, my $reverse_editing, my $supressPrint, my $revcom, my $border) = @_;
	$range = 1 unless $range; #default motif size = triplets
	my %triAllSeqs = ();
	my $triCount = 0;
	#Get editing sites from file
	my $sites_r; 
	if ($siteListFile =~ /clusters_/){ #get redundant sites from tabular cluster file
		(my $cluster_file, my $SorT) = split ('\|', $siteListFile); 
		$sites_r = redundantSiteList($cluster_file, $SorT); 
	}
	else{ #get non-redundant sites from siteList file
		$sites_r = sitesFromSiteList($siteListFile);
	}
	
	#Get triplets in editing sites of each seq in file
	my $inseq = Bio::SeqIO->new( -file => $seqFile, -format => 'fasta' );
	while ( my $seq = $inseq->next_seq ){
		#get triplets per seq
		my $s = lc $seq->seq(); 
		(my $coords) = $seq->id =~ /([^=]+:\d+-\d+[+-])/; 
		next unless exists $sites_r->{$coords}; #skip sequences that aren't in editing list.  
		my $triPerSeq = getTriplets($s, $acgt, $sites_r->{$coords}, $range, $reverse_editing); 
		
		#sum triplet frequencies 
		foreach my $t(keys %$triPerSeq){
			$triAllSeqs{$t} += $triPerSeq->{$t}; 
			$triCount += $triPerSeq->{$t}; #sum for normalization (for fraction; not used for enrichment)
		}
	}
	
	## normalize ##
	if ($normalize){ #norm by file or create fraction
		foreach my $t (keys %triAllSeqs){
			next if $triAllSeqs{$t} eq 'inf'; #***why should this happen?
			$triAllSeqs{$t} /= $triCount if $triCount;
			$triAllSeqs{$t} = nearest(.001, $triAllSeqs{$t}); #round to max of 3 digits after the decimal number
		}
		#get background triplet frequencies for normalization
		if (-e $normalize){ #seqfile was arged to use as background for enrichment calculation
			my $triBG; 
			my $fasta = $normalize; #change name for readability
			my $siteListArg = ($border ? $siteListFile : ""); 
			$triBG = getTripletsBG($fasta, $acgt, $range, 1, $siteListArg); 
			foreach my $t (keys %triAllSeqs){
				if ($triBG->{$t}){ #norm by background
					$triAllSeqs{$t} /= $triBG->{$t};
					$triAllSeqs{$t} = nearest(.001, $triAllSeqs{$t}); #round to max of 3 digits after the decimal number
				}
				else{ #avoid division by zero
					$triAllSeqs{$t} = "inf";
				}
			}
		}
		elsif($normalize ne '1'){ #normalize argument is filename that doesn't exist
			die "bad file name $normalize\n"; 
		}
	}
	#possibly revcom
	if ($revcom){
		my $revedHash = revcomHashKeys(\%triAllSeqs); 
		%triAllSeqs = %$revedHash; 
	}
	
	#print frequencies to file and stdin(optional) and return as hash
	my $tripletFile = $seqFile; 
	if (-e $normalize){ #background file for normalization was argumented
		$tripletFile =~ s/seqFasta/triplets/; 
	}
	elsif($normalize eq "1"){#fraction (relative frequency)
		$tripletFile =~ s/seqFasta/tripletsFrac/; 
	}
	else{ #frequency
		$tripletFile =~ s/seqFasta/tripletsFreq/; 
	}
	if ($border){
		$tripletFile =~ s/bordered_triplets/triplets/; 
	}
	$tripletFile =~ s/\.fa$/.txt/; 
	
	if ($range>1){
		my $filePrefix = "neighbors".$range; 
		$tripletFile =~ s/triplets/$filePrefix/; 
	}
	open (TRI, ">".$tripletFile) or die "open $tripletFile\n"; 
	foreach my $t (sort {$triAllSeqs{$b} <=> $triAllSeqs{$a}} keys %triAllSeqs){
		print $t ."\t". $triAllSeqs{$t} ."\n" unless($supressPrint); 
		print TRI $t ."\t". $triAllSeqs{$t} ."\n"; 
	}
	close(TRI); 
	
	return \%triAllSeqs; 
}

#Function: printTripletBackground2file - prints all triplets found in fasta file to freq file (count or frac)
#*** CHECK
#Input: 1. Fasta file (from which to get triplets)
#		2. acgt: center nuc of interest
#		3. siteList: the list of edited sites to border the regions from which to calculate frequency (from 1st edited to last edited + flanks of length 'range')
#		4. range: amount of sites before and after each potential editing site
#		5. normalize: return value format - 0 or "" (default): frequency; 
#										 1 = fraction of each triplet of all triplets 
#		6. tripletFile - output file
#		7. revcom: if to revcom the triplets
sub printTripletBG {
	(my $seqFile, my $acgt, my $siteListArg, my $range, my $normalize, my $tripletFile, my $revcom) = @_;
	my $triBG = getTripletsBG($seqFile, $acgt, $range, $normalize, $siteListArg);
	$triBG = revcomHashKeys($triBG) if($revcom); #revcom output (optional)
	if($tripletFile){
		open (TRI, ">".$tripletFile) or die "open $tripletFile\n"; 
	}
	else{
		open (TRI, \*STDIN) or die "open STDIN\n"; 
	}
	foreach my $t (sort {$triBG->{$b} <=> $triBG->{$a}} keys %$triBG){
		print TRI $t ."\t". $triBG->{$t} ."\n"; 
	}
	close(TRI); 
}

#Function: revcom all keys of a hash
sub revcomHashKeys {
	(my $hash_ref) = @_; 
	my $reved; 
	my %revedHash; 
	foreach my $key (keys %$hash_ref){
		$reved = revcom($key); 
		$revedHash{$reved} = $hash_ref->{$key}; 
	}
	return \%revedHash; 
}

#Function: revcom a string
#Note: all chars must be a,c,t,g (otherwise will be replaced with blank)
sub revcom {
	my $seq = shift; 
	my %map = ("a", "t", "c", "g", "g", "c", "t", "a", "A", "T", "C", "G", "G", "C", "T", "A"); 
	my @s = split ('',$seq); 
	for(my $i=0; $i<=$#s; $i++){
		$s[$i] = $map{$s[$i]} if exists $map{$s[$i]}; #special characters shouldn't be changed (such as Ns)
	}
	my $reved = join ('', @s); 
	return $reved;
}

#Function: get the frequency of each triplet with a specific nuc in center from file
#		Primary usage is for calculating background frequencies for motif inference of editing sites
#Input: 1. Fasta file (from which to calculate)
#		2. Range: amount of sites before and after each potential editing site
#		3. retfmt: return value format - 0 or "" (default): frequency; 
#										 1 = fraction of each triplet of all triplets 
#		4. siteList: the list of edited sites to border the regions from which to calculate frequency (from 1st edited to last edited + flanks of length 'range')
#RetVal: a hash with 1 level: HASH{triplet} = value (frequency, fraction or percent of each triplet (rounded to 3-digits after decimal))
sub getTripletsBG{
	(my $seqFile, my $acgt, my $range, my $retfmt, my $siteList) = @_;
	if (not -e $seqFile){ #file doesn't exist
		print "Error: File $seqFile doesn't exist\n"; 
		return 0; 
	}
	
	### get border per sequence, based on editing sites ### 
	my $borders_r = ""; 
	if ($siteList){ #get 0-base terminal editing sites as borders (see flank addition below)
		$borders_r = analysisSubs::getTerminalsFromSiteList($siteList);  
	}
	
	### get triplet frequencies for target nuc###
	my %triFreq = ();
	my $regex = "([acgt]{".$range."})" . $acgt . "([acgt]{".$range."})"; #create triplet (or quintuplet..) regex for parsing by regex
	my $inseq = Bio::SeqIO->new( -file => $seqFile, -format => 'fasta' );
	while ( my $seq = $inseq->next_seq ){
		my $s = lc $seq->seq();
		if ($borders_r){ #calc background from seq within borders of first and last editing sites (with immediate flanks of length $range)
			(my $coords) = $seq->id() =~ /([^=]+:\d+-\d+[+-])/;
			my $flank = $range; #seq used should include the flanks immediately next to first and last editing site. 
			my $start = $borders_r->{$coords}[0] - $flank; #include left flank of 1st editing site. 
			$start = 0 if $start < 0;
			my $end = $borders_r->{$coords}[1]; 
			my $substrLen = $end - $start + 1 + $flank; #substr length should include end site and right flank of last editing site. (only adds one flank because start is set to contain the other)
			$s = substr($s, $start, $substrLen); #(note: $start + $substrLen may be out-of-bounds - in that case substr returns only until the last char)
		}
		while ($s =~ /$regex/g){
			$triFreq{$&}++; 
			pos ($s) = pos ($s) - length($&) + 1; #move position back (start of last match + 1) so that regex doesn't skip overlapping triplets (e.g. "accc" is 2 triplets for "c" - "acc" and "ccc")
		}
	}
	#Alter output format - fraction instead of frequency (optional)
	if ($retfmt){
		my $sum = 0; 
		#sum triplet freqs
		foreach my $t (keys %triFreq){
			$sum += $triFreq{$t};
		}
		#convert from freqs to fractions
		foreach my $t (keys %triFreq){
			$triFreq{$t} /= $sum;
			$triFreq{$t} = nearest(.001, $triFreq{$t}); #round to max of 3 digits after the decimal number	
		}
	}
	return \%triFreq; 
}

#Function(helper): reverses the editing in a sequence
#Input: 1. sequence (lc/uc) 
#		2. site_ref: array_ref with editing positions 
#		3. pre_post: editing type to reverse (e.g. "GA" or "ga"; reverses A back to G) 
#		4. $retAr - flag (0=return seq as scalar, 1=return seq as array-ref).
#retVal: sequence as array-ref (default) or scalar with editing positions reversed
sub revEditing {
	(my $sequence, my $site_ref, my $pre_post, my $retAr) = @_;
	$sequence = lc $sequence;
	(my $pre_edit, my $post_edit) = split('', lc $pre_post);
	my @seq = split('',$sequence);
	foreach my $site (@$site_ref){
		if ($seq[$site-1] ne $post_edit){ #error: 'edited' nuc isn't the correct nuc
			print "Error: site should be \"$post_edit\" but isn't!\n"; 
			return 0;
		}
		$seq[$site-1] = $pre_edit; 
	}
	if ($retAr){ #return array-ref  
		return \@seq; 
	}
	else{ #Default: return sequence (scalar)
		return join ('', @seq); 
	}
}


#Function: parses sitelist file and returns a HoA with the site list ($sites{sequence ID} = array of sites) 
#Note: feature to add: id of sequence returned as coords only 
sub sitesFromSiteList {
	(my $siteListFile, my $fullDeflineKeys) = @_; 
	my %sites = (); 
	open(SITES, $siteListFile ) or die "siteList file $siteListFile didn't open\n";
	while (my $line = <SITES>)
	{ 
		(my $name)  = ( $line =~ /^(\S+)\t/ );
		(my @where) = ( $'    =~ /(\d+)/g );
		unless ($fullDeflineKeys){ #extract coords from name
			($name) = $name =~ /([^=]+:\d+-\d+[+-])/; 
		}
		$sites{$name} = \@where;
	}
	close(SITES);
	return \%sites; 
}

#Function: parses nuclist file and returns a HoA with the nuc list ($sites{sequence ID} = array of nucs in consensus) 
#Note: feature to add: id of sequence returned as coords only 
sub nucsFromNucsList {
	(my $nucListFile, my $fullDeflineKeys) = @_; 
	my %nucs = (); 
	open(NUCS, $nucListFile ) or die "nucList file $nucListFile didn't open\n";
	while (my $line = <NUCS>)
	{ 
		(my $name)  = ( $line =~ /^(\S+)\t/ );
		(my @whichNucs) = ( $'    =~ /([actg])/g );
		unless ($fullDeflineKeys){ #extract coords from name
			($name) = $name =~ /([^=]+:\d+-\d+[+-])/; 
		}
		$nucs{$name} = \@whichNucs;
	}
	close(NUCS);
	return \%nucs; 
}

#Function(helper): For each edited sequence get the positions of the terminal edited sites (first & last)
#Input: 1. siteList file.
#		2. base: 0/1 based sites (default 0-base)
#Retval: hash-of-array with borders(0/1 based, depending on input). 
#Note: 	Is used to calculate Nuc frequencies within edited regions.
sub getTerminalsFromSiteList {
	(my $siteList, my $base) = @_; 
	my %tSites = (); 
	my $sites_r = sitesFromSiteList($siteList); #get sites
	#calc terminals
	my @sites; 
	foreach my $coords (keys %$sites_r){
		@sites = @{$sites_r->{$coords}}; 
		my @terminals = ($sites[0], $sites[$#sites]); 
		unless ($base){
			$terminals[0]--;
			$terminals[1]--;
		}
		$tSites{$coords} = \@terminals; 
	}
	return \%tSites;
}



### Functions to find best parent ###
#Function: transform sequences in fasta file to 3-base
#Input: 1. input = filename or sequence (a,c,t,g,A,C,T,G,n,N only; will be converted to lc)
#note: converts to lower case
#Retrun: output file name
sub transformFileTo3base{
	(my $inFile, my $from_to, my $outFile) = @_;
	$outFile = $inFile .".". $from_to ."3base" unless $outFile; 
	(my $from, my $to) = split('', $from_to);
	my $inseq = Bio::SeqIO->new( -file => $inFile, -format => 'fasta' ); #infile
	my $outseq = Bio::SeqIO->new( -file => ">$outFile", -format => 'fasta' ); #outfile
	my $tseq; 
	while ( my $seq = $inseq->next_seq ){ 
		my $transedSeq = transformSeqTo3base($seq->seq(), $from_to); 
		$tseq = Bio::Seq->new( -seq => $transedSeq, -id  => $seq->id()); 
		$outseq->write_seq($tseq); 
	}
	return $outFile; 
}
#Function(helper): transform a sequence to 3-base
#Notes: 1. input isn't checked for compatibility 
#		2. converts to lc
sub transformSeqTo3base{
	(my $seq, my $from_to) = @_;
	$seq = lc $seq; 
	$from_to = lc $from_to; 
	(my $from, my $to) = split('', $from_to); 
	$seq =~ s/$from/$to/g;
	return $seq; 
}

#Function: find the best parent for each edited element. 
#         Other future options: "best parent" based on the parent element which minimizes the number of non G>A mismatches divided by alignment length: min(#non G>A mm / al_len)
#Input: 1. pairFile = coordinate editing-pairs (assumes that source and target are in last columns in file)
#		NOTE: A. if STreversed mode is on, then submit source and target files accordingly (e.g. if you want original GtoA to be AtoG then fastaS=A-seq and fastaT=G-seq)
#			  B. expected defline format is like in db (5 '=' separators)
#		2. fastaS = source file 
#		3. fastaT = target file
#		4. STreversed = target is one-before-last column and source is last (default is vice versa)
#		5. transform = if to transform to 3-base (e.g. GA (or ga) to transform G to A)
#Output: 
#		IMPORTANT: output retains order of source & target from original pair-file even if algorithm addresses them as reversed
#					However, if algorithm is reversed then prefix "STrev" will be added to output files
#		1. "bestPairs_" file = pairs of most similar squences in 3-base, selected from editing-pairs only
#		2. "otherBestPairs_" file = pairs that had equal bit-score to the pairs in bestPairs_ (where only the 1st was saved)
#Note: This is tricky in terms of Source/Target reversed. Be careful...
sub getBestSources {
	(my $pairFile, my $fastaS, my $fastaT, my $STreversed, my $transform, my $trimmed) = @_; 
	my $makeBestPairsClustersFile = 1; #CONST
	#create output filenames
	(my $path, my $filename) = $pairFile =~ /^(\S+\/)?([^\/]+)/; 
	my $suffix = $filename;
	$suffix =~ s/^[^_]+_// if ($filename =~ /_/); #leave only suffix (from first underscore; if no underscore leaves full name)
	my $rawClustersFile = $path . "clusters_" . $suffix; 
	$rawClustersFile =~ s/\.txt$/.tab/; 
	$suffix = ($STreversed ? "STrev_" : "") . $suffix; 
	$suffix = ($transform ? "3base_" : "")  . $suffix; 
	$suffix = ($trimmed ? "trimmed_" : "") . $suffix; 
	my $blastOut = $path . "t2sBLAST_". $suffix; #blast targets against sources
	my $bestPairsFile = $path . "bestPairs_" . $suffix;
	my $otherBestPairsFile = $path . "otherBestPairs_" . $suffix; 
	my $bestPairsClustersFile = $path . "bestPairsClusters_" . $suffix; 
	$bestPairsClustersFile =~ s/\.txt$/.tab/; 
	## BLAST ##
	#Transform seq files to 3-base (optional)
	my $tfastaS = $fastaS .".3base";
	my $tfastaT = $fastaT .".3base";
	if ($transform){
		transformFileTo3base($fastaS, $transform, $tfastaS); 
		transformFileTo3base($fastaT, $transform, $tfastaT); 
		$fastaS = $tfastaS; 
		$fastaT = $tfastaT; 
	}
	
	#Blastn to find best match in A-seq for each G-seq 
	my $log = $path . "makeblastdb.log.".$$.".txt"; #temp log file just to suppress output
	#system("makeblastdb -in $fastaS -dbtype nucl -logfile $log");
	#system("blastn -query $fastaT -db $fastaS -out $blastOut -strand plus -dust no -num_threads 8 -outfmt 6"); #-max_target_seqs 1
	system("formatdb -i $fastaS -p F -o T"); 
	system("blastall -p blastn -i $fastaT -d $fastaS -e 1e-50 -S 1 -F F -v 0 -a 8 -m 8 > $blastOut");	
	
	#del temp files
	if ($transform){ #delete temp 3-base files
		unlink ($tfastaS); 
		unlink ($tfastaT); 
	}
	unlink $log; #del makeblastdb logfile
	foreach my $suff (".nhr", ".nin", ".nsq"){ #del index files
		unlink ($fastaS . $suff); 
	}
	
	#Get all pairs from file and insert to hash_ref (format of hash: HASH{$target}{$source})
	my $t2s_r = getPairs($pairFile, $STreversed, 2); 
	#get best source for each target
	my %bestSource = ();
	my %otherBestSources = (); 
	my %bestScore = ();
	my $s; my $t; my $score; my $newT; 
	open (BLAST, $blastOut) or die "blastOut file $blastOut didn't open\n"; 
	while(my $l = <BLAST>){
		chomp $l; 
		my @fields = split (/\t/, $l);
		($t, $s, $score) = @fields[0,1,10]; #get target, source & bit-score
		($s) = $s =~ /([^=]+:\d+-\d+[+-])/; #get coords only, for comparison to pair file
		($t) = $t =~ /([^=]+:\d+-\d+[+-])/; #get coords only, for comparison to pair file
		if (exists $t2s_r->{$t}{$s}){ #pair was found as editing pair => candidate for best
			if (not exists $bestSource{$t}){ #1st result for the target - save it (BLAST results are sorted by bit-score for each target)
				$bestSource{$t} = $s;
				$bestScore{$t} = $score;
			}
			elsif($score == $bestScore{$t}){ #not 1st result for target but has best score too - save in separate file
				$otherBestSources{$t}{$s}=0; 
			}
		}
	}
	close(BLAST); 
	
	##print output (source\ttarget)
	#IMPORTANT: source and target here are as in original input file, despite the logic implemented in this function!
	#1. print best
	open (BEST, ">".$bestPairsFile) or die "couldn't open $bestPairsFile\n"; 
	foreach my $t (keys %bestSource){
		if ($STreversed){ #reverse output to retain source-target order from input pair-file
			print BEST  $t ."\t". $bestSource{$t} ."\n"; 
		}
		else{ #no need to reverse
			print BEST $bestSource{$t} ."\t". $t ."\n"; 
		}
	}
	close(BEST); 
	#2. print other best (tied in 1st place - save in separate file)
	open (OTHER, ">".$otherBestPairsFile) or die "couldn't open $otherBestPairsFile\n"; 
	foreach my $t (keys %otherBestSources){
		foreach my $s (keys %{$otherBestSources{$t}}){
			if ($STreversed){ #reverse output to retain source-target order from input pair-file
				print OTHER $t ."\t". $s ."\n"; 
			}
			else{ #no need to reverse
				print OTHER $s ."\t". $t ."\n"; 
			}
		}
	}
	close(OTHER); 
	
	#3. create clusters file for best pairs (for further analysis of only these clusters)
	if ($makeBestPairsClustersFile){
		open (my $c_raw_fh, $rawClustersFile) || die "open $rawClustersFile\n";
		open (my $c_best_fh, ">" . $bestPairsClustersFile) || die "open $bestPairsClustersFile\n"; 
		while(my $l = <$c_raw_fh>){
			chomp $l; 
			my @f = split(/\t/, $l); 
			my $s = ($STreversed ? $f[6] : $f[5]); 
			my $t = ($STreversed ? $f[5] : $f[6]); 
			if ($bestSource{$t} eq $s){
				print $c_best_fh $l ."\n";
			}
		}
		close($c_raw_fh); 
		close($c_best_fh); 
	}
	
	return \%bestSource;
}
#Function: blast3base - aligns one file to another as 3-base - output as table
#Input: qFa - query fasta file 
#		sFA - subject fasta file
#		blastOut - output file name
#		transform - what base to transform to which for 3-base (e.g. "ga" will transform g to a).
#		outfmt - output format of blast (default = 6 = tabular w.o. comments)
#Output: file 
sub blast3base {
	(my $qFa, my $sFa, my $blastOut, my $transform, my $outfmt, my $numThreads, my $maxTargets) = @_;
	#default parameters
	$transform = "ga" unless $transform; #transform g to a  
	$outfmt = 6 unless $outfmt; #tabular
	$numThreads = 10 unless $numThreads; #10 processors default
	
	#transform files to 3-base
	(my $path, my $filename) = $qFa =~ /^(\S+\/)?([^\/]+)/; 
	my $tqFa = transformFileTo3base($qFa, $transform); 
	my $tsFa = transformFileTo3base($sFa, $transform); 
	
	#Blastn to create tabular output
	my $log = $path . "makeblastdb.log.".$$.".txt"; #temp log file just to suppress output
	system("makeblastdb -in $tsFa -dbtype nucl -logfile $log");
	my $blastComm = "blastn -query $tqFa -db $tsFa -out $blastOut -strand plus -dust no -num_threads $numThreads -outfmt $outfmt"; 
	$blastComm .= " -max_target_seqs $maxTargets" if $maxTargets; 
	system($blastComm); #-max_target_seqs 1
	
	#del temp files
	unlink ($tqFa); 
	unlink ($tsFa); 
	unlink $log; #del makeblastdb logfile
	foreach my $suff (".nhr", ".nin", ".nsq"){ #del index files
		unlink ($tsFa . $suff); 
	}
	
}

#Function: blastAandGvsAluCons
#***under construction - for now just see in commands file
# sub blastAandGvsAluCons {
	# (my $qFa) = @_; 
	
	# my $qFa = $dir .""; 
	# my $sFa = "/home/alu/binknis/binknis_data/RepBase/ConsForMapping_18_8_13/SINE/Alu.fa"; 
	# my $blastOut = $dir . "seqA_blast3baseVsAllCons.txt";
	# analysisSubs::blast3base($qFa, $sFa, $blastOut); 

	# my $qFa2 = $qFa; 
	# $qFa2 =~ s/_A_/_G_/; 
	# $blastOut = $dir . "seqG_blast3baseVsAllCons.txt";
	# analysisSubs::blast3base($qFa2, $sFa, $blastOut); 
# }


#Function: gets pairs from pairs file and returns in hash
#Input: 1. pairFile = filename (must have last two columns of coords)
#		2. STreversed = if to reverse order of source target (0-no reverse (order is source\ttarget$), 1- reverse (order is target\tsource$)
#		3. retfmt = format of hashRef returned (0- pairs as keys, 1- source=>target, 2-target=>source)
#		4. retainFullName = if to return full names (1) or just coords (0) for pairs 
#Retval: ref to HoH: HASH{Target}{Source}
sub getPairs {
(my $pairFile, my $STreversed, my $retfmt, my $retainFullName) = @_;
	my %t2s = (); 
	my %s2t = (); 
	my %pairs = (); 
	my $source; 
	my $target; 
	#get pairs from file
	open(PAIRS, $pairFile) or die "$pairFile didn't open\n";
	while (my $line = <PAIRS>){
		chomp $line; 
		if ($STreversed){#reverse order - last col is source next-to-last is target
			($target, $source) = $line =~ /(\S+)\t(\S+)$/; 
		}
		else { #last col is source, next-to-last target
			($source, $target) = $line =~ /(\S+)\t(\S+)$/; 
		}
		unless($retainFullName){
			($source) = $source =~ /([^=]+:\d+-\d+[+-])/;
			($target) = $target =~ /([^=]+:\d+-\d+[+-])/;
		}
		
		$pairs{$source."\t".$target}=0; 
		$t2s{$target}{$source}=0; 
		$s2t{$source}{$target}=0; 
	}
	close(PAIRS); 
	
	return \%pairs if $retfmt==0;
	return \%s2t if $retfmt==1; 
	return \%t2s if $retfmt==2; 
}


### Mapping to consensus ###

##
#Function: getEditedNucsInCons - maps edited sites to the nucleotides found in consensus when mapped by blast
#***under construction 
#***The function of this subroutine is covered in getEditedPositionsInCons.
# sub getEditedNucsInCons {
	# (my $siteListFile, ) = @_; 
	 # my $sites_r = sitesFromSiteList($siteListFile);
# }


#Notes: 1. sitesListFile must be full path and present in Track-dir
sub getEditedPositionsInCons {
	#(my $siteListFile, my $idMapFile) = @_; 
	(my $siteListFile) = @_; 
	my $alwaysUseMostSimilarConsensus = 1; #CONST
	(my $dir, my $GA, my $suffix) = $siteListFile =~ /(\S+)\/siteList_([ACTG])_(\S+)\.txt/;
	my $subfamDir = $dir ."/SubfamFiles"; 
	(my $org, my $class, my $family, my $pval, my $th, my $control) = $suffix =~ /([^_]+)_([^_]+)_(\S+)_([^_]+)_([^_]+)(_control)?/; #needed only for class and family to access consensus file
	
	# my $rmskDir = "/home/alu/binknis/binknis_data/RepBase/ConsForMapping_18_8_13"; #CONST
	my $rmskRoot = "/home/alu/binknis/binknis_data/RepBase/ConsForMapping_12_1_14"; #CONST
	my $rmskDir = $rmskRoot; 
	### Create sequence files for each subfamily (G/A specific) and get subfamily names ###
	my $seqFile = "seqFasta_".$GA."_".$suffix.".fa";
	my $subfams = splitFastaBySubfam($dir, $seqFile, $GA);

	### get edited sites from file ###
	my $fullDeflineKeys = 1; #needed by getConsPosHist()
	my $sites_ref = sitesFromSiteList($siteListFile, $fullDeflineKeys);
	
	my $classConsFile = $rmskRoot ."/". $class .".fa"; 
	my $mostSimilarCons = getMostSimilarConsPerSubfam($subfamDir, $subfams, $GA, $classConsFile);
	
	### BLAST against the consensus sequence and create histogram of positions in consensus and list of nucleotides each element was mapped to ### 
	my %nucInConsList = (); 
	my $posConsFile = $dir ."/posCons_".$GA."_".$suffix.".txt";
	$rmskDir .= "/$class/$family";  
	open (POSCONS, ">$posConsFile") or die "$posConsFile didn't open\n";  
	my $unmappedFile = $subfamDir ."/". "unmapped_".$GA.".txt"; 
	my $mappedFile = $subfamDir ."/". "mapped_".$GA.".txt"; 
	open (UNMAPPED, ">$unmappedFile") or die "$unmappedFile didn't open\n";
	open (MAPPED, ">$mappedFile") or die "$mappedFile didn't open\n";
	foreach my $subfam (@$subfams){
		### run blast ###
		my $repbaseFile = $rmskDir ."/". $subfam.".fa"; #set direct file
		#if direct file doesn't exist - use most similar sequence based on blast search against class
		if (not -e $repbaseFile or $alwaysUseMostSimilarConsensus){ #construct repbase filename based on blast result
			(my $rmsk_sf, my $rmsk_fam, my $rmsk_class) = split('=', $mostSimilarCons->{$subfam}); 
			$repbaseFile = $rmskRoot ."/". join('/', $rmsk_class, $rmsk_fam, $rmsk_sf) .".fa"; 
		}
		#I commmented the option to use db_to_rmsk mapping for now
		# unless (-e $repbaseFile){ #if direct file doesn't exist then check if db_to_rmsk mapping contains a mapping for this subfam and reset file
			# my $mapping_file = "/home/alu/binknis/binknis_data/RepBase/files/mapping_output/name_mapping_files/db_to_rmsk_allClasses.txt"; 
			# open (my $map_fh, $mapping_file) or die "open $mapping_file\n"; 
			# while (my $l = <$map_fh>){
				# chomp $l; 
				# (my $db_sf, my $rmsk_sf) = split (/\t/, $l);
				# if ($db_sf eq $subfam){
					# $repbaseFile = $rmskDir ."/". $rmsk_sf.".fa"; 
					# last; 
				# }
			# }
			# close($map_fh); 
		# }
				
		if (-e $repbaseFile){
			my $editedFile = $subfamDir ."/seqFasta_".$GA."_".$subfam.".fa";  
			my $blastOutFile = $subfamDir ."/seqFasta_".$GA."_".$subfam.".blast.out";  
			unlink $blastOutFile if -e $blastOutFile; #erase old file because runBlast.pl concatenates
			#system("perl516 Tools/runBlast.pl $editedFile $repbaseFile $blastOutFile blastn 4 1e-10 1"); #cores=4, pval=1e-20, -S=1 (only + strand)
			#system("perl516 Tools/runBlast.pl $editedFile $repbaseFile $blastOutFile blastn 4 1e-2 1"); #cores=4, pval=1e-2, -S=1 (only + strand)
			my $blastPval = "1e-20"; my $cores = 4; my $blastStrand = 1;
			system("blastall -p blastn -i $editedFile -d $repbaseFile -e $blastPval -a $cores -S $blastStrand -F F -v 0 > $blastOutFile"); #formatdb was pre-run
			# system("blastn -query $editedFile -subject $repbaseFile -out $blastOutFile -strand plus -dust no -culling_limit 1"); #blast only plus strand
			
			### parse and print results ###
			# get consensus length (needed for parsing)
			my $consStream = Bio::SeqIO->new( -file => $repbaseFile, -format => 'fasta' );
			my $consSeq = $consStream->next_seq; 
			my $consLen = $consSeq->length();
			# get and print pos hist to family
			(my $subfamHist, $nucInConsList{$subfam} ) = getConsPosHist($blastOutFile, $consLen, $sites_ref, $subfam); #even though sites_ref contains multiple subfams the results will be for a single subfam because that's what the blast-file contains
			print POSCONS $subfam ."\t", join(' ',@$subfamHist) ,"\n";
			##print results to specific subfamily file
			my $histFile = $subfamDir."/consPosHist_".$GA."_".$subfam.".txt";
			writeSubfamHist($histFile, $subfamHist, $consSeq);
			#print log - where this subfam was mapped to
			(my $mappedTo) = $repbaseFile =~ /^\S+\/(\S+)\.fa$/; 
			print MAPPED join ("\t", $org, $class, $family, $subfam, $repbaseFile, $mappedTo)."\n"; 
		}
		else{
			print UNMAPPED join ("\t", $org, $class, $family, $subfam, $repbaseFile, "notMapped")."\n"; 
		}
	}
	close(POSCONS); 
	close(UNMAPPED); 
	close(MAPPED); 
	
	##write which nuc in consensus each edited site was mapped to
	my $nucListFile = $siteListFile; 
	$nucListFile =~ s/siteList/nucList/;
	open (NUCLIST, ">$nucListFile") or die "$nucListFile didn't open\n";
	foreach my $subfam(sort keys %nucInConsList){
		foreach my $seq(sort keys %{$nucInConsList{$subfam}}){
			print NUCLIST $seq , "\t",  join(' ',@{$nucInConsList{$subfam}{$seq}}) , "\n";
		}
	}
	close(NUCLIST); 
}

#*** create new files: consPosList (positions per element)



#Function(Helper): splitFastaBySubfam - 
#creates subdir in 'dir' for the output. 
#creates output files for each subfamily in the format of 'dir/subfamdir/seqFasta_[GA]_[subfam].fa'
sub splitFastaBySubfam{
	(my $dir, my $seqFile, my $GA) = @_; 
	my $subfamDir = $dir ."/SubfamFiles"; 
	mkpath $subfamDir;
	my %outStreams = ();
	my $inseq = Bio::SeqIO->new( -file => $dir ."/$seqFile",    -format => 'fasta' );
	while ( my $seq = $inseq->next_seq ){
		my $id = $seq->display_id(); 
		(my $subfam) = $id =~ /=([^=]+)$/; #get subfam
		
		if (not exists $outStreams{$subfam}){
			$outStreams{$subfam} = Bio::SeqIO->new( -file => ">".$subfamDir."/seqFasta_".$GA."_".$subfam.".fa",  -format => 'fasta' );
		}
		$outStreams{$subfam}->write_seq($seq);
	}
	my @subfams = sort keys %outStreams;
	return \@subfams; 
}

#Function: getMost similar consensus sequence per subfamily (needed when exact name of db subfam isn't in rmsk embl file)
# choose best consensus per subfam by choosing the subject most commonly producing best bitscore (first result).
sub getMostSimilarConsPerSubfam{
	(my $subfamDir, my $subfams, my $GA, my $classConsFile) = @_; 
	my %similarCons = (); 
	my %mostSimilarCons = ();
	#find candidates for best subfam per sequence
	foreach my $sf (@$subfams){
		#Blast subfam sequences to class's consensus sequences (parameters used: -F F disables filtering low complexity regions; -m tabular output format; -b 1 shows only best alignment;)
		#Only the best db sequence is shown for each query (-b 1); but will work without this filter too
		my $editedFile = $subfamDir ."/". "seqFasta_".$GA."_".$sf.".fa";
		my $blastOutFile = $subfamDir ."/"."blastToClass_". $GA ."_". $sf .".tab"; 
		my $blastPval = "1e-20"; my $cores = 4; my $blastStrand = 1;
		system("blastall -p blastn -i $editedFile -d $classConsFile -e $blastPval -a $cores -S $blastStrand -F F -v 0 -m 8 -b 1 > $blastOutFile"); 
		open (my $blast_fh, $blastOutFile) or die "open $blastOutFile\n"; 
		my $prev=""; 
		while(my $l = <$blast_fh>){
			chomp $l; 
			my @f = split(/\t/, $l);
			next if $f[0] eq $prev; #choose only first HSP for each sequence
			$prev = $f[0]; 
			my $sseqid = $f[1];
			$similarCons{$sf}{$sseqid}++; #save best hit
		}
		close($blast_fh); 
	}
	
	#choose best consensus per subfam by choosing the subject most commonly producing best bitscore (first result).
	foreach my $sf (keys %similarCons){
		my $maxCount=0; 
		foreach my $sseqid (keys %{$similarCons{$sf}}){
			if ($similarCons{$sf}{$sseqid} > $maxCount){
				$mostSimilarCons{$sf} = $sseqid;
				$maxCount = $similarCons{$sf}{$sseqid};
			}
		}
	}
	return \%mostSimilarCons;
}

#***possibly need to add a "next" command if is the incorrect subfamily
#Notes: 1. Alignments to gapped are skipped, not counted. 
#		2. All sequences must be aligned to the same sequence (all mapping will be inserted into a single hash)
#		3. sites_ref must have full defline keys, not only coords, for subfam to be extracted (if $subfam flag is off, it isn't necessary)
#		4. IDs in output will be identical to those in sitesList file. 
#Note: The original algorithm was flawed, but fixed to avoid double mapping in different HSPs. 
#***	Would be a good idea to double-check that the index array indMapped works precisely.

sub getConsPosHist{
	(my $blastFile, my $consLen, my $sites_ref, my $subfam) = @_; 
	my @posHist = (0) x $consLen;
	my %sites_ref_coords = (); 
	foreach my $defline (keys %$sites_ref){
		(my $coords) = $defline =~ /([^=]+:\d+-\d+[+-])/; 
		$sites_ref_coords{$coords} = $sites_ref->{$defline}; 
	}
	##create datastructures for saving nucs in consensus and avoiding double mapping of same position twice
	my %consNuc = (); #to which nucleotide each editing site was mapped to in consensus
	my %numToMap = (); #num of sites left to map
	foreach my $seq (keys %$sites_ref){
		next if ($subfam and not subfamFromID($seq, $subfam)); #subfam flag is on and seq isn't from the subfam
		my @arr = (0) x scalar(@{$sites_ref->{$seq}});
		$consNuc{$seq} = \@arr;
		$numToMap{$seq} = scalar(@arr);
	}
	##map to consensus
	my $in = new Bio::SearchIO( -format => 'blast', -file => $blastFile );
	while ( my $result = $in->next_result )
	{
		my $seqID = $result->query_name();
		next if ($subfam and not subfamFromID($seqID, $subfam)); #subfam flag is on and seq isn't from the subfam
		(my $coords) = $seqID =~ /([^=]+:\d+-\d+[+-])/; 
		unless(exists $consNuc{$seqID}){ #Find correct format that was in sitesList file (if not full defline)
			my $seqID_no_number = $seqID;
			$seqID_no_number =~ s/^\d+=//;
			if (exists $consNuc{$seqID_no_number}){ #Format produced from runClusterFinder after changing default output to tabular (Dec. 2013)
				$seqID = $seqID_no_number; 
			}
			elsif(exists $consNuc{$coords ."=". $subfam}){ #Format of minimal data for this subroutine. Hasn't been used yet. 
				$seqID = $coords ."=". $subfam; 
			}
			else{
				# print Dumper(\%consNuc); #*** 			
				# print "seqID: " . $seqID ."\n"; #***
				# print "seqID no number: " . $seqID_no_number ."\n"; #***
				die "Sequence deflines extracted from sites_ref don't match any expected format - exiting getConsPosHist for $blastFile. SeqID example: $seqID\n"; 
			}
		}
		# print $subfam ."\t". $coords ."\n"; 
		my @indices = @{ $sites_ref_coords{$coords} };
		my @indMapped = (0) x $#indices; #keep track of indices that were mapped
		#print $result->query_name(), ": @indices \n";
		while ( my $hit = $result->next_hit )
		{ 
			while ( my $hsp = $hit->next_hsp )
			{
				my $strq      = $hsp->query_string;
				my $strs      = $hsp->hit_string;
				my $alignment = $hsp->homology_string;
				my $q_start   = $hsp->start('query'); #1-base
				my $q_end   = $hsp->end('query'); 
				my $s_start   = $hsp->start('subject');
				my @sq        = split( //, $strq );
				my @ss        = split( //, $strs );
				my @al        = split( //, $alignment );
				my $al_len = $hsp->length('total'); 
				my $q_gaps    = 0;
				my $s_gaps    = 0;
				my @q2al = (0) x $al_len;
				my @al2s = (0) x $al_len;
				#map query sequence pos to alignment pos and map alignment pos to subject sequence pos 
				foreach my $i ( 0 .. $#al )
				{
					if ( $sq[$i] eq '-' )#count query gaps
					{
						$q_gaps++;
						next;
					}
					if ( $ss[$i] eq '-' )#count subject gaps
					{
						$s_gaps++;
						next;
					}
					#map query str to pos in align and then pos in align to subject str
					$q2al[$q_start + $i - $q_gaps - 1] = $i; #maps 0-base in q to 0-base in al
					$al2s[$i] = $s_start + $i - $s_gaps - 1; #maps 0-base in al to 0-base in s
				}
				foreach my $i (0 .. $#indices)
				{	
					next if $indMapped[$i]; #skip - was already mapped
					my $ind = $indices[$i] - 1; #input sites are 1-base so decrement
					if ($indices[$i] >= $q_start and $indices[$i] <= $q_end) #skip indices that aren't in hsp range
					{
						if ($ss[$q2al[$ind]] ne '-'){ #skip if aligned to gap
							#print "q_ind = $ind, al_ind = $q2al[$ind], s_ind = $al2s[$q2al[$ind]]\n";					
							$posHist[$al2s[$q2al[$ind]]]++; #increment pos in histogram
							$consNuc{$seqID}->[$i]= $ss[$q2al[$ind]]; #map this index in q to the nuc it aligned to in s
							$numToMap{$seqID}--; #one nuc mapped - one less to map for this seq
							$indMapped[$i]=1; #this index was mapped - don't map it again
						}
						else{ 
							# print "gap\n"; 
						}
					}
					last unless $numToMap{$seqID};
				}
				last unless $numToMap{$seqID};
			}
			last unless $numToMap{$seqID}; 
		}
	}
	return (\@posHist, \%consNuc); 
}

#Function: subfamFromID - check if a sequence is from a specific subfamily by it's id (defline)
#Input: id - format where family is last in line after a "=" 
#		subfam - the subfam to check if it belongs to
#Retval: if found (0/1)
sub subfamFromID{
	(my $id, my $subfam) = @_; 
	my $found = 0; 
	(my $id_sf) = $id =~ /\S+=(\S+)/; 
	$found = 1 if $id_sf eq $subfam; #subfam found as is
	$id_sf =~ s/\//_/;
	$id_sf =~ s/\?//; 
	$found = 1 if $id_sf eq $subfam; #subfam found with modified format (no question marks and underscore instead of slashes)
	return $found; 
}


#Function(Helper): writeSubfamHist - changes the original vector presentation of edited sites (posCons file) to a histogram of editing presentation ()
sub writeSubfamHist{
	(my $histFile ,my $subfamHist, my $consSeq) = @_; 
	open( HIST, ">".$histFile ) or die "didn't open $histFile\n";
	my @str = split( //, $consSeq->seq() );
	my %mapping = ('a'=>1,'c'=>2,'g'=>3,'t'=>4);
	foreach my $i ( 0 .. $#str )
	{
		print HIST $i + 1, "\t"; #index of pos in consensus
		#amount of sites mapped to this pos
		if ( $subfamHist->[$i+1] > 0 ) 
		{
			print HIST "$subfamHist->[$i+1]\t";
		}
		else
		{
			print HIST "0\t";
		}
		#the code of nucleotide in the pos
		if ($str[$i] =~ /[acgt]/) 
		{
			print HIST $mapping{$str[$i]} , "\n";
		}
		else
		{
			print HIST "5\n";
		}
	}
	close(HIST);
}

#Function: nucListToFreq - converts a nucList file to count and fraction of nucleotide that each nucleotide was mapped to
sub nucListToFreq {
	(my $nucListFile) = @_; 
	my $nucFreqFile = $nucListFile; 
	$nucFreqFile =~ s/nucList/nucListFreq/;
	my %map = ('a'=>0,'c'=>1,'g'=>2,'t'=>3);
	
	#get count and fractions
	my $nuc_ref = nucsFromNucsList($nucListFile, 1); 
	my %nucFreq = (); 
	my %total = (); 
	foreach my $seq (keys %$nuc_ref){
		my @initAr = (0) x 4; 
		$nucFreq{$seq} = \@initAr;
		foreach my $n (@{$nuc_ref->{$seq}}){
			$nucFreq{$seq}[$map{$n}]++; 
			$total{$seq}++; 
		}
	}	
	
	#print output
	open (NUCFREQ, ">". $nucFreqFile) or die "open $nucFreqFile\n"; 
	foreach my $seq (keys %$nuc_ref){
		print NUCFREQ $seq;
		foreach my $count (@{$nucFreq{$seq}}){
			print NUCFREQ "\t" . $count; 
		}
		foreach my $count (@{$nucFreq{$seq}}){
			my $frac = nearest(.01, $count / $total{$seq}) if $total{$seq}; 
			print NUCFREQ "\t" . $frac; 
		}
		print NUCFREQ "\n";
	}
	close (NUCFREQ); 
}

#Function: create a Freq for source/target sequence for each pair 
sub nucFreqPerPairs {
	(my $clustFile, my $SorT, my $GA) = @_; 
	$GA = "G" unless $GA; 
	
	#Map sites to consNuc for all sites in source/target
	(my $path, my $suffix) = $clustFile =~ /^(\S+)clusters_(\S+)\.tab$/; 
	my $siteListFile = $path ."siteList_".$GA."_".$suffix.".txt"; 
	my $sl = sitesFromSiteList($siteListFile, 0);
	my $nucListFile =  $path ."nucList_".$GA."_".$suffix.".txt";
	my $nl = nucsFromNucsList($nucListFile, 0);
	
	#Get sites per pair
	my $nucListPerPairFile =  $path ."nucListPerPair_".$GA."_".$suffix.".txt";
	my %pn = (); #pair nucs
	my $ps = getPairSites($clustFile, $SorT);
	
	#get nucList for sites per pair only
	foreach my $pair (sort keys %$ps){
		(my $s, my $t) = split(/\t/, $pair); 
		my $coords = ($SorT ? $t : $s); 
		my $i=0; #pair-sites index
		my $j=0; #all seq-sites index
		my @p_sites = @{$ps->{$pair}};
		my @p_nucs = (); 
		# print "p_sites: ", join(" ", @p_sites), "all_sites: ", join(" ", @{$sl->{$coords}}),"\n"; 
		while($i <= $#p_sites){
			while ($sl->{$coords}[$j] != $p_sites[$i]){
				# print $sl->{$coords}[$j] ." ". $p_sites[$i] ."\n";
				if ($j == scalar(@{$sl->{$coords}})) { die "error for $clustFile, passed array limits in nucFreqPerPairs\n";} #just for debugging, shouldn't happen
				$j++; 
			}
			push(@p_nucs, $nl->{$coords}[$j]); 
			$i++;
		}
		$pn{$pair} = \@p_nucs; 
	}
	
	open (NLPPF, ">". $nucListPerPairFile) or die "open $nucListPerPairFile\n"; 
	foreach my $pair (sort keys %pn){
		my $pairOut = $pair; 
		$pairOut =~ s/\t/|/;
		print NLPPF $pairOut ,"\t", join(' ',@{$pn{$pair}}), "\n"; 
	}
	close(NLPPF); 
	#Create freq file per pairs
	nucListToFreq($nucListPerPairFile);
}

#Function: getRedundantSiteList - same as sitesFromSiteList
#SorT - sites in source or in target (source is defined as the first ) (0- source, 1- target)
sub redundantSiteList {
	(my $tabClustFile, my $SorT) = @_; 
	##Get unique edited sites between each pair
	my $pairsSites = getPairSites($tabClustFile, $SorT); 
		
	##Create list of unique for each source(/target)
	my %sitesPerSeq = (); #redundant-sites per source(/target)
	foreach my $p (keys %$pairsSites){
		(my $source, my $target) = split(/\t/, $p); 
		my $SorTseq = ($SorT ? $target : $source); #sequence to use for hash returned (default: source)
		push( @{$sitesPerSeq{$SorTseq}}, @{$pairsSites->{$p}} ); 
	}
	return \%sitesPerSeq; 
}

#Function(Helper): getPairSites - gets a hash of all editing sites between a pair of sequences
#if cluster file has multiple entries of same pair the unique sites in their alignment will be returned. 
#Helper for markEditingSitesInBLAST
# SorT = 0 - get for source, 1 - get for target
#Fields for cluster file input parsing: [0]mmType, [1]assembly, [2]class, [3]family, [4]subfam, [5]coordsS, [6]coordsT, [7]num_edited_sites, [8]num_all_mms, [9]total_prob, [10]whereS, [11]whereT, [12]num_clusts, [13]mmSerials, [14]clusters_span_woGaps, [15]clusters_span, [16]alignment_len, [17]mmCount_str
sub getPairSites {
	(my $tabClustFile, my $SorT) = @_; 
	my %ps = (); #pair-sites
	open (TC, $tabClustFile) or die "open $tabClustFile\n"; 
	while(my $l = <TC>){
		next unless $l =~ /\S+/;
		chomp $l;
		my @f = split(/\t/, $l);
		my $s = $f[5];
		my @s_sites; 
		unless ($SorT){ #source
			@s_sites = split(',', $f[10]);
		}
		else{ #target
			@s_sites = split(',', $f[11]);
		}
		my $t = $f[6];
		my $pair = $s ."\t". $t; 
		push (@{$ps{$pair}}, @s_sites); 
	}
	close(TC);
	
	foreach my $p ( keys %ps){ #get unique sites per editing-pair (not per edited element)
		my @tmp = sort {$a <=> $b} uniq(@{$ps{$p}});
		$ps{$p} = \@tmp; 
	}
	
	return \%ps; 
}

#Function: getNucsForPairs - get hash of nucleotides found in a position relative to each editing site in cluster-pairs
#Input: 
#	1. tabClustFile - full path to tabular cluster file
#	2. SorT - get nucs in source (0) or target (1); Source and Target nucs are defined in the subroutine
#	3. 
#	4. pairsFilterFile - file with pairs to retain (needed when working in a track directory which is after filtering)
#						This file needs to be graph_ file from the Track directory.
#***Modify for mismatches other than GA and CT [see (? : )]
# Notes: 
#	1. For filters you must specify the tabClustFile after filtering
sub getNucsForPairs {
	(my $tabClustFile, my $SorT, my $position, my $pairsFilterFile) = @_; 
	my $sNuc = "G"; 
	my $tNuc = "A"; 
	my %pairNucs = ();
	##get sites per pair
	#get pairs if clusters should be filtered for a set of pairs (if pairsFilterFile with pairs to retain was specified)
	my $pairsFromFilter; 
	if ($pairsFilterFile){
		$pairsFromFilter = getPairs($pairsFilterFile);
	}
	my $pair_sites = getPairSites( $tabClustFile, $SorT); 
	my %seq2ps = (); #hash{seq(source or target)}{pair}{sites-per-pair}
	foreach my $p (keys %$pair_sites){
		if ($pairsFilterFile){ #filter - if file with pairs to retain was specified
			next unless exists $pairsFromFilter->{$p}; 
		}
		my @seqs = split(/\t/, $p); 
		$seq2ps{$seqs[$SorT]}{$p} = $pair_sites->{$p}; #SorT is used here as a flag if to set keys as source or target
	}
	
	##get nucs relative to sites of pairs
	#create seq file name to fetch nucs from AND create outfile name
	(my $path, my $suffix) = $tabClustFile =~ /^(\S+\/)?clusters(_[^\/]+)\.tab$/;
	if ($pairsFilterFile){ 
		($path, $suffix) = $pairsFilterFile =~ /^(\S+\/)?graph(_[^\/]+)\.txt$/;
	}
	my $seqFile = $path ."seqFasta_". ($SorT ? $tNuc : $sNuc) . $suffix .".fa";
	my $outFile = $path ."nucListPairs_" .($SorT ? $tNuc : $sNuc) ."_pos". $position. $suffix .".txt";
	#fetch nucs from seqFile
	open (OUT, ">".$outFile) || die "open $outFile\n"; 
	my $inseq = Bio::SeqIO->new( -file => $seqFile, -format => 'fasta' );
	while ( my $seq = $inseq->next_seq ){
		(my $id) = $seq->id =~ /([^=]+:\d+-\d+[+-])/;
		if (exists $seq2ps{$id}){ 
			# print $id .":\n"; #***
			foreach my $pair ( keys %{$seq2ps{$id}} ){
				# print $pair ."\t"; #***
				$pairNucs{$pair} = (); 
				foreach my $site ( @{$seq2ps{$id}{$pair}} ){
					my $wantedNuc = substr($seq->seq, $site-1+$position, 1); 
					push ( @{$pairNucs{$pair}},  $wantedNuc); 
					# print $wantedNuc ." "; #***
				}
				#print output per pair: 
				my $outPair = $pair; 
				$outPair =~ s/\t/,/; 
				# ***my $tmpLen = scalar(@{$pairNucs{$pair}}); 
				#***print OUT $pair ,"\t", "len: ",$tmpLen, "\n";
				print OUT $outPair ,"\t", join (" ", @{$pairNucs{$pair}}), "\n";
			}
		}
	}
	close(OUT); 
	return (\%pairNucs, $outFile); 
}

#Function(helper): uniq - Return sorted & unique elements in array 
sub uniq {
    return sort keys %{{ map { $_ => 1 } @_ }};
}

#Function: createConsensus - creates a consensus sequence from the sequences in a fasta file using MSA
#*** To do 
sub createConsensus{
	(my $fastaFile) = @_; 
	#***use some package that creates MSA and creates consensus from it.
}

#Function: - makeStarredAlign
#*** To do 
sub markEditingSitesInBLAST{
	(my $blastFile, my $clusterFile) = @_; 
	my $in = new Bio::SearchIO(-format => 'blast', -file   => $blastFile);
	
	# while 

}

### Subroutine alignment_to_mm_arr ###
#Function: parses an hsp and returns arrays of all the mismatches in the alignment 
#Input: a searchIO hsp object
#Returns: 3 HoHoA refs: 	1. mm_q: 12 arrays containing all positions in query of every mismatch. 
#							2. mm_s: same, for subject. 
#							3. ind: 12 arrays containing serial numbers of each type of mismatch ('serial number' = position in alignment)
sub alignment2mmArray{
	(my $hsp_ref) = @_; 
	my $hsp = $$hsp_ref;
	# my @nucs = qw/a c g t/; 
	# print lc $hsp->query_string , "\n"; #***
	my @q = split( //, lc $hsp->query_string );
	my @s = split( //, lc $hsp->hit_string );
	my $alignment = $hsp->homology_string;
	
	#position hashes for query and subject sequences
	my %mm_q = (); 
	my %mm_s = (); 
	my %ind = (); 
	my %q2gaps = (); 
	my %s2gaps = (); 
	# for my $nuc1 (@nucs){
		# for my $nuc2 (@nucs){
			# next if $nuc1 eq $nuc2; 
			# $mm_q{$nuc1}{$nuc2} = (); 
			# $mm_s{$nuc1}{$nuc2} = (); 
			# $ind{$nuc1}{$nuc2} = (); 
		# }
	# }
	#position and counter init
	my $q_start  = $hsp->start('query');
	my $s_start  = $hsp->start('subject');
	my $mm_count = 0;
	my $q_gaps   = 0;
	my $s_gaps   = 0;
	
	while ( $alignment =~ /( )/g )
	{
		my $p = ( pos $alignment ) - 1; #pos function is 1-base, so needs to be decremented to fit $s and $q which are 0-base. 
		if ( $q[$p] eq '-' ){ #gap in query
			push(@{$s2gaps{$s[$p]}}, $s_start + $p - $s_gaps); #save nuc and position in subject that aligned to gap in query
			$q_gaps++;
		}
		elsif ( $s[$p] eq '-' ){ #gap in sbjct
			push(@{$q2gaps{$q[$p]}}, $q_start + $p - $q_gaps); #save nuc and position in query that aligned to gap in subject
			$s_gaps++;
		}
		else{ #mismatch
			push (@{ $mm_q{$q[$p]}{$s[$p]} } , $q_start + $p - $q_gaps); 
			push (@{ $mm_s{$s[$p]}{$q[$p]} } , $s_start + $p - $s_gaps); 
			$mm_count++;
			push (@{ $ind{$q[$p]}{$s[$p]} } , $mm_count); 
		}
	}	
	return (\%mm_q, \%mm_s, \%ind, \%q2gaps, \%s2gaps);
}

#Function: getBlastMMfreqs - prints the distribution of each type of mismatch in all alingments in a BLAST file (or directory)
#To read from STDIN, specify 'STDIN' as input
#Input: 
#	$blastFile - input file
#	$outFile - output file
#	$coordsOut - flag (0- print defline as is, 1-extract coords) for "1", must have coords at head of line or after '=' or '|' chars 
#	$append - flag if to append to output (1) file or override (0) 

#Output fields: 
#[0]SourceName [1]TargetName [2]$hsp_count,$hsp->hsp_length,$hsp->evalue,$hsp->score,$totalMMs [3]ga|ct|gc|gt|ca|ta|ag|tc|cg|tg|ac|at 
#Explanation: spaces are tabs. 

sub calcBlastMMfreqs{
	(my $blastFile, my $outFile, my $coordsOut, my $append) = @_; 
	print $blastFile ."\n";  #***
	my @nucs = qw/a c g t/;
	my @mmOrder = qw/ga ct gc gt ca ta/;
	
	#set input file
	my $gzipped = ($blastFile =~ /\.gz$/ ? 1 : 0); 
	my $in;
	if ($blastFile eq 'STDIN'){
		$in = new Bio::SearchIO(-format => 'blast', -fh   => \*STDIN);
	}
	elsif($gzipped){
		open(my $pipeFromZip, "gunzip -c $blastFile |") or die "open gzipped $blastFile\n"; 
		$in = new Bio::SearchIO(-format => 'blast', -fh   => $pipeFromZip);
	}
	else{
		$in = new Bio::SearchIO(-format => 'blast', -file   => $blastFile);
	}
	
	#set output file handle
	my $out_fh; 
	if(not $outFile){ #defualt - use input file as prefix for output and add ".mmFreqs.txt"
		$outFile = $blastFile .".mmFreqs.txt" unless ($outFile); 
		open ($out_fh, ">" . ($append ? ">" : "") . $outFile) or die "open $outFile\n";
	}
	elsif ($outFile eq 'STDOUT'){  #write to stdout
		$out_fh = \*STDOUT; 
	}
	else{ #write to specified file
		open ($out_fh, ">" . ($append ? ">" : "") . $outFile) or die "open $outFile\n";
	}
	
	#read blast file
	while (my $result = $in->next_result){
		while(my $hit = $result->next_hit){
			my $hsp_count=0; 
			next if $result->query_name() eq $hit->name(); #skip alignments between sequence to itself
			#set query and subject names for output
			my $outNames; 
			if ($coordsOut){
				(my $qCoords) = $result->query_name() =~ /([^=|]+:\d+-\d+[+-])/;
				(my $sCoords) = $hit->name() =~ /([^=|]+:\d+-\d+[+-])/;
				$outNames = join("\t", $qCoords , $sCoords); 
			}
			else{
				$outNames = join("\t", $result->query_name() , $hit->name()); 
			}
			#get and write mismatch array and alignment stats to output
			while(my $hsp = $hit->next_hsp){
				$hsp_count++;
				# print "HSP number $hsp_count\n"; #***
				#parse alignment to mismatch arrays
				(my $mm_q_r, my $mm_s_r, my $ind_r, my $q2gaps_r, my $s2gaps_r)  = alignment2mmArray(\$hsp);
				# print scalar(@{$mm_q_r->{"g"}{"a"}}) ."\t". $mm_s_r ."\t". $ind_r ."\n" if $mm_q_r->{"g"}{"a"}; #***
				#print output
				my @mmCount = (0) x 12; #hold number of each type of mismatch (as in @mmOrder, and then the 6 respective reverse mismatches; e.g. GA and then AG)
				my $nuc1; my $nuc2; 
				foreach my $i (0 .. 5){ #six mismatches, as in @mmOrder
					# print "mmorder: "; print "@mmOrder\n"; #***
					($nuc1, $nuc2) = split('', $mmOrder[$i]);
					# print $nuc1 ."\n"; 
					# print "hello $nuc1 $nuc2" .  $mmOrder[$i] ."\n"; #***
					$mmCount[$i] = scalar(@{$mm_q_r->{$nuc1}{$nuc2}}) if ($mm_q_r->{$nuc1}{$nuc2}); #mismatch count
					# print "scalar\n" . scalar(@{$mm_q_r->{$nuc1}{$nuc2}}) . "\n" if ($mm_q_r->{$nuc1}{$nuc2}); #mismatch count#***
					$mmCount[$i + 6] = scalar(@{$mm_q_r->{$nuc2}{$nuc1}}) if ($mm_q_r->{$nuc2}{$nuc1}); #reverse mismatch count
					# print "scalar\n" . scalar(@{$mm_q_r->{$nuc2}{$nuc1}}) . "\n"  if ($mm_q_r->{$nuc2}{$nuc1}); #reverse mismatch count#***
				}
				my $totalMMs = 0; 
				# print "1: "; 
				# print "@mmCount\n";
				$totalMMs += $_ for @mmCount; 
				# print "2: "; 
				# print "@mmCount\n";
				
				##Save nucleotide composition of query and subject (in aligning region)
				#Get nuc composition of query string (region aligned)
				my $nc_ref_q = getSeqNucComposition($hsp->query_string()); 
				my @ncCountQ = (0) x 4; 
				foreach my $i(0 .. 3){
					$ncCountQ[$i] = $nc_ref_q->{$nucs[$i]} if $nc_ref_q->{$nucs[$i]}; 
				}
				
				#Get nuc composition of query string (region aligned)
				my $nc_ref_s = getSeqNucComposition($hsp->hit_string()); 
				my @ncCountS = (0) x 4; 
				foreach my $i(0 .. 3){
					$ncCountS[$i] = $nc_ref_s->{$nucs[$i]} if $nc_ref_s->{$nucs[$i]}; 
				}
								
				##Save count of number of times each nuc aligned to a gap
				#number of times each query nuc aligned to a gap
				my @gapCountQ = (0) x 4; 
				foreach my $i(0 .. 3){ #query aligned to gap count per nucleotide
					$gapCountQ[$i] = scalar(@{$q2gaps_r->{$nucs[$i]}}) if ($q2gaps_r->{$nucs[$i]}); 
				}
				
				#number of times each subject nuc aligned to a gap
				my @gapCountS = (0) x 4; 
				foreach my $i(0 .. 3){  #subject aligned to gap count per nucleotide
					$gapCountS[$i] = scalar(@{$s2gaps_r->{$nucs[$i]}}) if ($s2gaps_r->{$nucs[$i]}); #gap count
				}
				
				#Build output
				my $alignData = join (',', $hsp_count, $hsp->hsp_length, $hsp->evalue, $hsp->score, $totalMMs); 
				my $mmData = join ('|', @mmCount); 
				my $gapStrQ = join('|', @gapCountQ); 
				my $gapStrS = join('|', @gapCountS); 
				my $nucStrQ = join('|', @ncCountQ);
				my $nucStrS = join('|', @ncCountS);
				
				print $out_fh join("\t", $outNames, $alignData, $mmData, $gapStrQ, $gapStrS, $nucStrQ, $nucStrS) , "\n";
			}
		}
	}
	
	if ($outFile ne 'STDOUT'){
		close($out_fh); 
	}
}

# wrapper for getBlastMMfreqs - for whole directory
#$out : 0 - directory (will create file for each subfam) or 1 - one file
#Note: added $dataDir
sub calcBMMForg{
	(my $dataDir, my $org, my $class, my $out) = @_;
	my $home = $ENV{"HOME"};
	my $dirPrefix = "BlastStats"; 
	unless ($out){ #default - create output directory in home directory 
		$out = $home . "/". $dirPrefix; 
		mkdir $out; 
	}
	my $coordsOut = 0; 
	my $multiFile = (-d $out ? 1 : 0); #if output file (0) or directory (1) specified in input
	my $append = ($multiFile ? 0 : 1); #append only if it's not multifile
	my $blastDir = join ('/', $dataDir , $org, $class, "results", "blasts"); 
	my $fams = getAll::families($dataDir, $org, $class); 
	foreach my $fam (@$fams){
		my $outFamDir = join ('/', $out, $dirPrefix, $org, $class, $fam); 
		mkpath $outFamDir if $multiFile; 
		my $subfamFiles = getAll::name_files($dataDir, $org, $class, $fam);
		foreach my $sff (@$subfamFiles){
			my $blastFile = join ('/', $blastDir, $fam, $sff .".gz");
			my $outFile; 
			$outFile = ($multiFile ? $outFamDir ."/". $sff : $out); 
			calcBlastMMfreqs($blastFile, $outFile, $coordsOut, $append); 
		}
	}
}

#*** Need to check if this works *** 
#$mm: two letters resembling the mismatch that we'd like to check its motif (e.g. GA is from G to A); lc/uc are fine. 
#$name: full defline or regex that matches the defline
#Note: runClusterFinder enables
sub motif_per_alignment{
	(my $blastFile, my $name1, my $name2, my $mm) = @_; 
	my @nucs = qw/a c g t/;
	my %nucMap = ('a'=> 0, 'c' => 1, 'g' => 2, 't' => 3); 
	(my $from, my $to) = split('', lc $mm); #extract nuc 'from' and 'to' from mismatch arg. 
	my $in; 
	if ($blastFile eq 'STDIN'){
		$in = new Bio::SearchIO(-format => 'blast', -fh   => \*STDIN);
	}
	else{
		$in = new Bio::SearchIO(-format => 'blast', -file   => $blastFile);
	}
	#read blast file
	while (my $result = $in->next_result){
		while(my $hit = $result->next_hit){
			if ($result->query_name() =~ /$name1/ and $hit->name() =~ /$name2/){
				(my $org, my $class, my $fam, my $subfam) = $result->query_name =~ /^\d+=([^=]+)=[^=]+=([^=]+)=([^=]+)=([^=]+)$/; 
				my $hsp_count=0; 
				while(my $hsp = $hit->next_hsp){
					$hsp_count++; 
					# print "HSP number $hsp_count\n"; 
					#parse alignment to mismatch arrays
					(my $mm_q_r, my $mm_s_r, my $ind_r)  = alignment_to_mm_arr(\$hsp); 
					if ($mm_q_r->{$from}{$to}){
						#get query string and erase gaps
						my $q_seq = $hsp->query_string; 
						$q_seq =~ s/-//g;
						print $q_seq ."\n"; 
						my $start = $hsp->start('query'); 
						my @mm_arr = @{$mm_q_r->{$from}{$to}}; 
						for(my $i=0; $i<=$#mm_arr; $i++){
							$mm_arr[$i] = $mm_arr[$i] - $start + 1;
						}
						#get nuc frequency/motif and print output
						my $range = 2; 
						my $freq_ref = getNucFrequencyPerPos($q_seq, \@mm_arr, $range, 0, 0);
						my $fam_nuc_freq = getAll::nucFreqSubfam($dataDir, $org, $class, $fam, $subfam); 
						# my @poses = (-$range .. -1, 1 .. $range);
						my @poses = (-$range .. $range);
						foreach my $pos (@poses){
							print $pos ."\t"; 
							my %sum_freq_ref = (); #for final normalization (after normalization by nucleotide frequency, we want to normalize so that each positions sums to 100%)
							foreach my $nuc (@nucs){ #normalize by frequency of nucleotides in subfamily
								$freq_ref->{$pos}{$nuc} /= $fam_nuc_freq->[$nucMap{$nuc}];
								$sum_freq_ref{$pos} += $freq_ref->{$pos}{$nuc}; 
							}
							foreach my $nuc (@nucs){ 
								$freq_ref->{$pos}{$nuc} /= $sum_freq_ref{$pos} unless $freq_ref->{$pos}{$nuc} == 0; #*** Later try to understand if this normalization is inadequate. i.e. is it an enrichment score? if so it shouldn't be normalized to sum-1.
								print $freq_ref->{$pos}{$nuc} * 100 . ($nuc eq 't' ? "\n" : "\t"); 
							}
						}
					}
				}
			}
		}
	}
}

#Function: getSeqNucComposition - Returns the nucleotide composition of a single sequence
#Input: 1. Sequence (will be converted to lc)
#		2. retfmt: return value format - 0 or "" (default): frequency; 
#										 1 = fraction of each nuc in seq
#		4. startEnd: The start and end positions between which to calc the nucComposition (0-base; seperated with pipe '|').
#Retval: nc - a,c,g,t,g+a,c+t,c+g,a+t frequencies or fractions
sub getSeqNucComposition {
	(my $seq, my $retfmt, my $startEnd) = @_;
	my @nucs = ('a', 'c', 'g', 't'); 
	my @nucPairs = ('ga', 'ct', 'cg', 'at'); 
	my %nc = ();
	foreach my $n (@nucs, @nucPairs){
		$nc{$n}=0;
	}
	
	$seq = lc $seq;
	if ($startEnd){ #get substr
		(my $start, my $end) = split(',', $startEnd); 
		$seq = substr($seq, $start, ($end - $start + 1)); 
	}
	$seq =~ s/[^actg]//g;
	my $sum = length($seq); 
	#calc compistion
	$nc{$_}++ for (split("", $seq)); #nuc
	foreach my $np (@nucPairs){ #nuc pairs
		(my $first, my $second) =  split('',$np); 
		$nc{$np} = $nc{$first} + $nc{$second};
	}
	
	if ($retfmt){
		#sum nuc freqs in each pos
		#convert from freqs to fractions
		foreach my $n (@nucs, @nucPairs){
			$nc{$n} /= $sum unless $sum==0; #norm and avoid division by 0
			$nc{$n} = nearest(.001, $nc{$n}); #round to max of 3 digits after the decimal number	
		}	
	}
	return \%nc;
}

#Function: getAllSeqNucComposition - returns the nucleotide composition of all sequences in a fasta file
sub getAllSeqNucComposition {
	(my $seqFile, my $retfmt, my $siteList) = @_;
	my @nucs = ('a', 'c', 'g', 't');
	my @nucPairs = ('ga', 'ct', 'cg', 'at');
	my %ncPerSeq = (); #per sequence
	my %ncAll = (); #all sequences summed
	foreach my $nuc (@nucs, @nucPairs){
		$ncAll{$nuc}=0;
	}
	## get border per sequence, based on editing sites
	my $borders_r = ""; 
	if ($siteList){
		$borders_r = analysisSubs::getTerminalsFromSiteList($siteList);
	}
	
	### fetch nucleotide frequencies around each nucleotide ###
	my $inseq = Bio::SeqIO->new( -file => $seqFile, -format => 'fasta' );
	while ( my $seqObj = $inseq->next_seq ){
		my $seq = lc $seqObj->seq();
		
		my $coords; 
		my $noCoords; 
		if ($seqObj->id() =~ /([^=|]+:\d+-\d+[+-])/){ #coords found in defline
			($coords) = $seqObj->id() =~ /([^=|]+:\d+-\d+[+-])/;
			$noCoords=0; 
		}
		else{ #no coords present in defline - can't use borders
			$coords = $seqObj->id();
			$noCoords=1; 
		}
		my $ncOneSeq; 
		if ($borders_r and not $noCoords){ #get substr within first and last editing site -- in order to calc background from seq within borders of first and last editing sites
			my $start = $borders_r->{$coords}[0];
			$start = 0 if $start < 0;
			my $end = $borders_r->{$coords}[1];
			($ncOneSeq) = getSeqNucComposition($seq, 0, ($start.','.$end)); 
		}
		else{
			($ncOneSeq) = getSeqNucComposition($seq, 0); 
		}
		
		#Get and sum frequencies per sequence
		for my $n (@nucs, @nucPairs){
			$ncAll{$n} += $ncOneSeq->{$n}; 
			$ncPerSeq{$coords}{$n} = $ncOneSeq->{$n}; 
		}
	}
	
	### Normalize 
	if($retfmt){
		#norm for all-seqs
		my $sum = 0; 
		$sum += $ncAll{$_} for (@nucs); 
		foreach my $n (@nucs, @nucPairs){
			$ncAll{$n} /= $sum unless $sum==0; #norm and avoid division by 0
			$ncAll{$n} = nearest(.001, $ncAll{$n}); #round to max of 3 digits after the decimal number	
		}	
		
		#norm for per-seq
		foreach my $s(keys %ncPerSeq){
			$sum = 0;
			$sum += $ncPerSeq{$s}{$_} for (@nucs); 
			foreach my $n (@nucs, @nucPairs){
				$ncPerSeq{$s}{$n} /= $sum unless $sum==0; #norm and avoid division by 0
				$ncPerSeq{$s}{$n} = nearest(.001, $ncPerSeq{$s}{$n}); #round to max of 3 digits after the decimal number
			}	
		}
		
	}
	
	return (\%ncAll, \%ncPerSeq);
}

#Function that prints the sequence nucleotide composition to file
sub printAllSeqNucComposition{
(my $seqFile, my $retfmt, my $siteListFile, my $outFile) = @_; 
	my @nucs = ('a', 'c', 'g', 't');
	my @nucPairs = ('ga', 'ct', 'cg', 'at');
	(my $ncAll, my $ncPerSeq) = getAllSeqNucComposition($seqFile, $retfmt, $siteListFile);
	
	#print nuc comp for all seqs together
	open (OUT, ">".$outFile) or die "open $outFile\n"; 
	my $outStr = "ALL";
	foreach my $n(@nucs, @nucPairs){
		$outStr .= "\t" . $ncAll->{$n};
	}
	print OUT $outStr."\n"; 
	
	#print nuc comp for each seq
	foreach my $coords (sort keys $ncPerSeq){
		$outStr = $coords; 
		foreach my $n (@nucs, @nucPairs){
			$outStr .= "\t" . $ncPerSeq->{$coords}{$n};
		}
		print OUT $outStr."\n"; 
	}
	close(OUT);
}


#Function: given a sequence and list of modified sites, returns the modifications made in each site.
#	$mm = 2 letter mismatch (e.g. GA = G-to-A)
#	$pos_r = positions mutated and the number of occurences: $pos_r->{$pos} = Num-changes
#			Poses are 0-base
#	$orfStart & $orfEnd =  start and end positions of ORF within the sequence (optional; default - whole sequence is the ORF).
#			start/end are 0-base.
#	$codons_r = (optional) - send codon list for efficiency.
#	To do (Optionally) -1.  add parameter: $isPost = sequence is post-editing (by default, sequence is pre-editing)
#	Note: does not deal with exons/introns. So ORF needs to be uninterrupted by gaps.
sub getCodonChangesPerSeq {
	(my $mm, my $seq, my $pos_r, my $orfStart, my $orfEnd, my $codons_r) = @_; 
	#convert to uc for hash to work
	$mm = uc $mm; 
	(my $pre_mm, my $post_mm) = split(//, $mm); 
	$seq = uc $seq; 
	my @s = split(//, $seq); #seq to array of 1 nt in each cell
	
	$orfStart = 0 unless $orfStart; #default start = 0
	$orfEnd = $#s unless $orfEnd; #default end = last pos in seq (0-base)

	#Set codon hash.
	my %cods; 
	if ($codons_r){ #if sent as parameter for efficiency
		%cods = %$codons_r; 
	}
	else{
		%cods = (
		'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A', 'TTA' => 'L',
		'TTG' => 'L', 'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L',
		'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R', 'AGA' => 'R',
		'AGG' => 'R', 'AAA' => 'K', 'AAG' => 'K', 'AAT' => 'N', 'AAC' => 'N',
		'ATG' => 'M', 'GAT' => 'D', 'GAC' => 'D', 'TTT' => 'F', 'TTC' => 'F',
		'TGT' => 'C', 'TGC' => 'C', 'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P',
		'CCG' => 'P', 'CAA' => 'Q', 'CAG' => 'Q', 'TCT' => 'S', 'TCC' => 'S',
		'TCA' => 'S', 'TCG' => 'S', 'AGT' => 'S', 'AGC' => 'S', 'GAA' => 'E',
		'GAG' => 'E', 'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
		'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G', 'TGG' => 'W',
		'CAT' => 'H', 'CAC' => 'H', 'TAT' => 'Y', 'TAC' => 'Y', 'ATT' => 'I',
		'ATC' => 'I', 'ATA' => 'I', 'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V',
		'GTG' => 'V', 'TAA' => 'X', 'TAG' => 'X', 'TGA' => 'X'); #X represents "Stop".
	}
	
	##RetVal data structures
	my %changeSum = (); #Syn/NS/StopGain/StopLoss/StartLoss/Inconsistent.
	my %posSum = (); #Syn/NS/StopGain/StopLoss/StartLoss/Inconsistent.
	my @changeTypes = ("Syn", "NS", "StopGain", "StopLoss", "StartLoss", "Inconsistent"); 
	for my $c (@changeTypes){
		$changeSum{$c} = 0; 
		$posSum{$c} = 0; 
	}
	
	#Get changes
	foreach my $p (keys %{$pos_r}){
		### skip positions outside of ORF
		next unless $p >= $orfStart and $p <= $orfEnd; 
		##1. get codon pre and post editing 
		my $cpos = ($p - $orfStart + 1) % 3; #codon pos (0/1/2)
		my $pre; 
		my $post; 
		if($s[$p] ne $pre_mm){ #inconsistent (the position isn't the expected pre-editing nt)
			$changeSum{"Inconsistent"} += $pos_r->{$p}; 
			$posSum{"Inconsistent"}++; 
			next; 
		}
		else{ #Consistent (modified position contains the expected nucleotide
			if($cpos == 1){ #1st-pos in codon
				$pre = join('', $pre_mm, $s[$p+1], $s[$p+2]);
				$post = join('', $post_mm, $s[$p+1], $s[$p+2]);
			}
			elsif($cpos == 2){ #2nd-pos in codon
				$pre = join('', $s[$p-1], $pre_mm, $s[$p+1]);
				$post = join('', $s[$p-1], $post_mm, $s[$p+1]);
			}
			else{ #3rd-pos in codon
				$pre = join('', $s[$p-2], $s[$p-1], $pre_mm);
				$post = join('', $s[$p-2], $s[$p-1], $post_mm);
			}
		}
		
		#2. get type of mismatch
		if(not exists $cods{$pre}){ #Problematic character in sequence (Maybe N)
			next; 
		}
		my $mm_type; 
		if($cods{$pre} ne $cods{$post}){ #Non-Synonymous
			if($cods{$post} eq 'X'){ #premature stop codon
				$mm_type = 'StopGain'; 
			}
			elsif($cods{$pre} eq 'X'){#stop loss
				$mm_type = 'StopLoss'; 
			}
			elsif($p < $orfStart + 3 and $cods{$pre} eq 'M'){ #start-codon loss
				$mm_type = 'StartLoss'; 
			}
			else{ #Non-synonymous
				$mm_type = 'NS'; 
			}
		}
		else{ #Synonyomous
			$mm_type = "Syn"; 
		}
		#save mismatch type mutation-sum and position-sum
		$changeSum{$mm_type} += $pos_r->{$p}; 
		$posSum{$mm_type}++; 
	}
	return(\%changeSum, \%posSum);
}

#Test getCodonChangesPerSeq subroutine
sub testCodonChangesPerSeq{
	use Data::Dumper;
	my $seq = "ATATGATGAAACCTCCTTCTCAATGGTAG"; #ATG AAA CCT CCT TCT CAA TGG TAG
	my %posFreq = ( '7' => 15, '5' => 1, '19' => 2, '25' => 7, '21' => 3, '22' => 2); 
	# print Dumper(\%posFreq); exit; 
	(my $changeSum_r, my $posSum_r) = getCodonChangesPerSeq("GA", $seq, \%posFreq, 5, length($seq)); 
	print "****Changes sum****\n"; 
	print Dumper($changeSum_r);
	print "****Positions sum****\n"; 
	print Dumper($posSum_r);
}

#Function: Get 
#Output: orfStats file
#Note: 	1. skips proteins that have introns (or for some other reason are interrupted (otherwise, the downstream function getCodonChangesPerSeq can generate errors)
#		2. IMPORTANT: Not currently built to run generically... (only consts, no input args)
sub getConsAndPosHists {
	my $mm = "GA"; 
	my $GA = "A"; 
	my $filter = "lowConfidence_byClusters_AmoreDiv_mostGMapGinPairs_bestSources"; 
	# my $suffix = "Zebrafinch_LTR_ERVL_1e-0_5"; #$organism . "_" . $class . "_" . $family . "_" . $pval . "_" . $th
	my $suffix = "Zebrafinch_LTR_ERVK_1e-0_5"; #$organism . "_" . $class . "_" . $family . "_" . $pval . "_" . $th
	
	my $dir = "/home/alu/binknis/Data/Zebrafinch/LTR/results/Tracks/Filtered/tracks_".$suffix."/".$mm."/".$filter;
	my $posConsFile = $dir ."/". "posCons_".$GA."_".$suffix.".txt"; 
 	my $cdsPosFile = "/home/alu/binknis/binknis_data/RepBase/RepBase18.02.embl.withExtractedProteins/proteins_10_9_14/All.Repbase.ref.proteins_cdsPositions.txt"; 
	my $mapFile = $dir ."/". "SubfamFiles/mapped_".$GA.".txt";
	my $orfStatsFile = $dir ."/". "orfStats_".$GA."_".$suffix.".txt"; 
	
	my %poses = ();
	my %cdss = (); 
	my %sf2cons = ();
	my %consFiles = (); 
	
	##Get Positions per subfam Hash
	open (POS, $posConsFile) or die "open $posConsFile\n"; 
	while(my $l = <POS>){
		chomp $l; 
		(my $subfam, my @posArr) = split(/\s+/, $l);
		for(my $i=0; $i<=$#posArr; $i++){
			if ($posArr[$i] > 0){
				$poses{$subfam}{$i} = $posArr[$i]; 
			}
		}
	}
	close(POS);
	
	##Get ORF mapping
	open (CDS, $cdsPosFile) or die "open $cdsPosFile\n"; 
	while(my $l = <CDS>){
		chomp $l; 
		(my $pname, my $cons, my $species, my $pstart, my $pend, my $exonNum) = split(/\t/, $l); 
		next if $exonNum != 0; #SKIP PROTEINS WITH INTRONS (the code is 0-no introns, >0 - gives exon serial number for those with exons).
		$cdss{$cons}{$pname} = $pstart .",". $pend; 
	}
	close(CDS); 
	
	##Get ORF mapping
	open (MAPPING, $mapFile) or die "open $mapFile\n"; 
	while(my $l = <MAPPING>){
		chomp $l; 
		(my $org, my $class, my $fam, my $subfam, my $consFile, my $cons) = split(/\t/, $l);
		$sf2cons{$subfam} = $cons; #subfam -> consensus to which it mapped
		$consFiles{$cons} = $consFile; #consensuses -> respective sequence files
	}
	close(MAPPING);
	
	##Get ORF stats and print output
	(my $org, my $class, my $fam, my $pval, my $th) = split('_', $suffix); 
	open(ORFOUT, ">".$orfStatsFile) or die "open $orfStatsFile\n"; 
	foreach my $subfam (sort keys %poses){
		my $cons = $sf2cons{$subfam}; 
		foreach my $prot (sort keys %{$cdss{$cons}}){
			#get ORF's DNA sequence
			my $inseq = Bio::SeqIO->new( -file => $consFiles{$cons}, -format => 'fasta' );
			my $seqObj = $inseq->next_seq;
			my $consSeq = uc $seqObj->seq();
			(my $orfStart1base, my $orfEnd1base) = split(',', $cdss{$cons}{$prot});
			(my $changeSum_r, my $posSum_r) = getCodonChangesPerSeq($mm, $consSeq, $poses{$subfam}, $orfStart1base - 1, $orfEnd1base - 1);
			foreach my $mmType (sort keys %$changeSum_r){
				my $outstr = join("\t", $org, $class, $fam, $subfam, $cons, $prot, $mmType, $changeSum_r->{$mmType}, $posSum_r->{$mmType}) . "\n"; 
				print ORFOUT $outstr; #print output
			}
		}
	}
	close(ORFOUT);
	
}


#Function: creates a fasta file from a list of requested sequences from the DB (in Data dir)
#			The script is a basic step in createTrackFiles2.pl and used to be a script getFastaFromCoords.pl. 
#Input: coordinate file (can be multiple format, but must have coords as 'chr:start-end' and afterwards 'subfamily')
#		column numbering is 0 based 
# $fastaOut : can provide with or without suffix (.fa will be the suffix; and _subfam.fa if file per subfam is activated)
#			  Can contain full path!
# $outFormatFlag : 0 - full, 1 - only coords, 2 - "coords=subfam"
# $filePerSubFamOutput : 0 - one large file. 1 - a file per subfamily
# $columnsToGet : specifies the column containing the coords; 0-base.
# $coordsToSubfam : HASH - subfams to use based on coords (needed for case where DB contains 1 file per fam). 0 - DB contains file per subfam, so no need to change subfam names. 
#Output: Fasta file with sequences (1 per family if flag is on)
#Notes: 1. coord filename must contain the suffix created by createTrackFiles.pl (Org, class, fam, pval, th)
#		2. subfam names with slashes and/or question-marks will have them retained in deflines but are substituted with underscores in filenames. 
#		3. Added dataDir as arg (instead of ../Data)
# use Data::Dumper; 
#*** ADD option to read from zipped file

sub getFastaFromCoords{
	(my $dataDir, my $coordFile, my $fastaOut, my $columnToGet, my $outFormatFlag, my $filePerSubFamOutput, my $coordsToSubfam) = @_; 
	unless(-d $dataDir){
		print "Error: invalid dir arg to analysisSubs::getFastaFromCoords. First arg must be dataDir!\n";
	}
	
	# print Dumper($coordsToSubfam);
	my $coords; 
	my $coordToGet; 
	my $org, my $class, my $fam; 
	my %subfamToCoords = (); 
	#extract Organism, Class and family from file name
	if ($coordFile =~ /\S+[\/]*[^_\/]+_([^_\/]+)_([^_\/]+)_(\S+)_(1e-\d+)_(\d+)(_control)?.txt$/){
		($org, $class, $fam) = ($coordFile =~ /\S+[\/]*[^_\/]+_([^_\/]+)_([^_\/]+)_(\S+)_(1e-\d+)_(\d+)(_control)?.txt$/); 
	}
	else{
		die "insufficient data in file name: $coordFile - need Organism, Class and Family\n"; 
	}	

	#create hash with the subfamily and coords wanted ({subfam}{coords})
	open (COORDS, "$coordFile") or die "$coordFile didn't open\n"; 
	while (my $line = <COORDS>){
		chomp $line; 
		# print "line: " . $line ."\n"; #*** 
		my @fields = split (/\s+/,$line); 
		$coordToGet = $fields[$columnToGet];
		# print "coord to get: " . $coordToGet ."\n"; #***
		($coords, my $subfam) = $coordToGet =~ /([^=]+:\d+-\d+[+-])\S*=(\S+)$/;
		# print "coords - $coords , subfam - $subfam "."\n"; #***
		if($coordsToSubfam){ #replace subfam found in defline with subfam in $coordsToSubfam HASH
			$subfam = $coordsToSubfam->{$coords}; 
		}
		$subfam =~ s/\?//g;  
		$subfam =~ s/\//_/g;
		$subfamToCoords{$subfam}{$coords} = 0;
	}
	close (COORDS);

	#copy the wanted sequences from the subfamily files in the DB
	my $outStream;
	$fastaOut =~ s/\.fa$//; #to avoid adding .fa to existent .fa suffix
	$outStream = Bio::SeqIO->new( -file => ">".$fastaOut.".fa",  -format => 'fasta' ) unless $filePerSubFamOutput;
	foreach my $subfam (sort keys %subfamToCoords){
		my $dbFile = $dataDir."/".$org."/".$class."/db/files_".$fam."/Seq_".$subfam; 
		if($coordsToSubfam){ #only one file per family in input DB
			$dbFile = $dataDir."/".$org."/".$class."/db/files_".$fam."/Seq_".$fam; 
		}
		if ($filePerSubFamOutput){
			$outStream = Bio::SeqIO->new( -file => ">".$fastaOut."_".$subfam .".fa",  -format => 'fasta' ); #subfamily-specific output file
		}
		copySubfamSequences($dbFile, $outStream, $subfamToCoords{$subfam}, $outFormatFlag, $coordsToSubfam); 
	}
}

#Function: reads the fasta db file and copies only those existing in the coords hash to the output stream
#Is helper function for 
#	The label is modified with respect to the label(id) flag
sub copySubfamSequences {
	(my $dbFile, my $faOut, my $coords_ref, my $idFlag, my $c2sf) = @_;
	my $inseq = Bio::SeqIO->new( -file => $dbFile,    -format => 'fasta' );
	while ( my $seq = $inseq->next_seq ){
		my $id = $seq->display_id(); 
		(my $coords) = $id =~ /^\d+=[^=]+=(\S+:\d+-\d+[+-])/; #get coords from id
		if (exists $coords_ref->{$coords}){ #write the sequence to output if the sequence's coords exist in the input file
			if($c2sf){ #replace subfam from defline with subfam from coordsToDefline hash
				$id =~ s/[^=]+$//; 
				$id .= $c2sf->{$coords}; 
				$seq->display_id($id);
			}
			#modify label, unless flag == 0
			if ($idFlag == 1){ #label of "coords=subfam"
				(my $subfam) = $id =~ /\S*=(\S+)$/;
				$seq->display_id($coords ."=". $subfam); 
			}
			elsif ($idFlag == 2){ #label of "coords"
				$seq->display_id($coords); 
			}
			#write seq to output
			$faOut->write_seq($seq);
		}
	}
}


return 1; 