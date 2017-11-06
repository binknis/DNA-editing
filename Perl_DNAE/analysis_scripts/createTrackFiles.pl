#Function: Parses a specific cluster file and creates several UCSC-track and analysis-output files
#Input: Specification of a cluster file (see args)
#	subdir = in result file
#	filter = select only part of sequences in cluster file: Format := [pairs/seqs/seqsG/seqsA/subfam/subfams]|[1/0]|[bed-file/gff-file/pair-file/subfam/subfam-file]|description. example: "seqsA|1|bed|this_describes_filter" (means: a sequence file for A sequences | retain only these | BED format input | "this_describes_filter" will be appended to output dir)
#     Filter explanation and pairing: 
	#		1. pairs - 3rd param will be two columns in "chrN:start-end[+-]" format 
	#		2. seqs - filter for both G and A seqs (BED or GFF file)
	#		3. seqsG - filter for G seqs 
	#		4. seqsA - filter for A seqs 
	#		5. subfam - one subfam
	#		6. subfams - a file with a list of subfams		
	#   1/0: retain (1) or remove (0) the specified sequences
	# Note: Input format of sequences in filter-file must be full(with multiple "=" signs) or simply coords (with or without strand)
#		Filters already completed: Maybe all... need to verify ***
#	*** IMPORTANT: MAKE SURE FILTER APPLIES TO ALL SCRIPTS RUN FROM THIS SCRIPT *** I think I already did so...
#	#
#
#
#
#
#Output: creates - 	1. Tabular cluster file; 
#					2. graph file (0-base starts); 
#					3. graph2 - fields: org class fam name assembly source_coords target_coords (Note: assumes name of source and target are the same)
#					3. track files (gff files are 1-base start); 
#					4. 
#					5. 
#exec-example:  perl analysis_scripts/createTrackFiles.pl Human LINE L1 1e-6 10  0 Clusters_80_percent_length_G &
use strict;
use lib $ENV{"HOME"} ."/Perl_DNAE";
use getAll; 
use lib $ENV{"HOME"} ."/Perl_DNAE/analysis_scripts"; 
use analysisSubs; 


( my $organism, my $class, my $family, my $pval, my $th, my $control, my $subdir, my $filter) = @ARGV;
$pval = "1e-" . $pval if $pval =~ /^\d+$/;
$subdir = 0 if $subdir eq ''; 
$control = 0 if $control eq '';
my $pmotif = '1e-3';
my $filter_subdir = ($filter ? "/Filtered" : ""); 
my $f_type;  my $f_negate; my $f_pairs=""; my $f_seqsG=""; my $f_seqsA=""; my $f_subfams=""; #filter variables used later
#### build file name and track dir ####
#file name 

my $cluster_file = "../Data/" . $organism . "/" . $class . "/results"; 
$cluster_file .= "/$subdir" unless ($subdir eq 0 or $subdir eq ""); 
$cluster_file .= "/clusters_" . $organism . "_" . $class . "_" . $family . "_" . $pval . "_" . $th;
#track dir
my $dir = "../Data/" . $organism . "/" . $class . "/results"; 
$dir .= "/$subdir" unless ($subdir eq 0 or $subdir eq ""); 
$dir .= "/Tracks"; 
$dir .= $filter_subdir;
mkdir $dir; #make "Tracks" dir; later will create specific 'tracks' dir
$dir .= "/tracks_" . $organism . "_" . $class . "_" . $family . "_" . $pval . "_" . $th;
#file and dir trailers
if ($control)
{
	$cluster_file .= "_control";
	$dir     .= "_control";
}
if ($filter){
	(my $f_suffix) = $filter =~ /([^\|]+)$/; 
	$dir .= "_".$f_suffix;
	$filter =~ s/\|[^\|]+$//; 
}

(my $suffix) = $cluster_file =~ /clusters_(\S+)$/;
my $tabular_file = $dir ."/clusters_". $suffix .".tab"; 
$cluster_file .= ".txt";
#don't parse if cluster's file is empty (avoids creating empty track files). 
exit if (-z $cluster_file); 
mkdir $dir;
if ($filter){
	createFilterFile($dir, $filter); #Create filter-description file
	($f_negate, $f_pairs, $f_seqsG, $f_seqsA, $f_subfams) = parseFilter($filter);
}

## Table variables ##
my $assembly; my $class; my $family; my $name;
my $chrG; my $strandG; my $startG; my $endG;
my $chrA; my $startA; my $endA; my $strandA;
my $mmSerials; my $locG=""; my $locA="";
my $align_length; my $direct_mms; my $all_mms; my $prob; my $num_clusts; 

## DSs for graph and tracks ## 
my $coordsG; 
my $coordsA; 
my %coords_G_all  = ();
my %coords_A_all  = (); 
my %minProb_G = (); 
my %minProb_A = (); 
my %sites_G = (); 
my %sites_A = (); 
my %graphEntries = ();
my %coordsToDefline = (); 

#### Parse cluster file and create tabular cluster file ###
my $retained = 0; #for filter count (to exit of 0 are retained after filter)
open (CLUSTS ,"<" . $cluster_file) || die ("couldn't open $cluster_file"); 
open (TAB ,">" . $tabular_file) || die ("couldn't open $tabular_file"); 
my $tuple;
while(my $line  = <CLUSTS>)
{
	if($line =~ /^Found/)
	{
		### READ CLUSTER ### 
		my $g_name_l  = <CLUSTS>; chomp $g_name_l; 
		my $a_name_l = <CLUSTS>; chomp $a_name_l; 
		my $prob_clustNum_l = <CLUSTS>;	chomp $prob_clustNum_l; 
		my $mm_serial_num_l = <CLUSTS>;	chomp $mm_serial_num_l; 
		my $locG_l = <CLUSTS>; chomp $locG_l; 
		my $locA_l = <CLUSTS>;	chomp $locA_l; 
		my $len_mmNum_l = <CLUSTS>;	chomp $len_mmNum_l; 
	
		### G name line ###
		#G name = 1082=hg19=chr1:144836998-144837209-=LINE=L1=HAL1
		($assembly, $chrG, $startG, $endG, $strandG, $class, $family, $name) = $g_name_l =~ /G name = \d+=([^=]+)=([^:]+):(\d+)-(\d+)([+-])=([^=]+)=([^=]+)=(\S+)/; 
		
		### A name line ###
		# A name = 118=hg19=chr1:16882278-16882479+=LINE=L1=HAL1
		($assembly, $chrA, $startA, $endA, $strandA, $class, $family, $name) = $a_name_l =~ /A name = \d+=([^=]+)=([^:]+):(\d+)-(\d+)([+-])=([^=]+)=([^=]+)=(\S+)/;
		
		### total_prob line ###
		# Total probability = 0.000172152888096637, 1 cluster
		($prob, $num_clusts) = $prob_clustNum_l =~ /Total probability = (\S+), (\d+) cluster/; 		
		
		### Edited MM serial line ###
		# Edited MM serial no.    :     1     2     3     4     5
		($mmSerials) = $mm_serial_num_l =~ /Edited MM serial no.    :\s+(.+\d+)/; 
		$mmSerials =~ s/\s+/,/g;
		
		### G locations line ###
		# Locations G seq:    72    77    89    99   141
		($locG) = $locG_l =~ /Locations G seq:\s+(.+\d+)/; 
		$locG =~ s/\s+/,/g;
		
		### A locations line ### 
		# Locations A seq:    71    76    88    98   140
		($locA) = $locA_l =~ /Locations A seq:\s+(.+\d+)/; 
		$locA =~ s/\s+/,/g;
		
		### Length, Mismatches ### 
		# Total length: 201, Direct mismatches: 5, All mismatches: 16
		($align_length, $direct_mms, $all_mms) = $len_mmNum_l =~ /Total length: (\d+), Direct mismatches: (\d+), All mismatches: (\d+)/;
		
		my $tuple = $assembly."\t".$class."\t".$family."\t".$name."\t"
					.$chrG."\t".$startG."\t".$endG."\t".$strandG."\t".$locG."\t"
					.$chrA."\t".$startA."\t".$endA."\t".$strandA."\t".$locA."\t"
					.$align_length."\t".$direct_mms."\t".$all_mms."\t".$mmSerials."\t".$prob."\t".$num_clusts; 
		
		$coordsG = $chrG .":". $startG ."-". $endG .$strandG; 
		$coordsA = $chrA .":". $startA ."-". $endA .$strandA; 
		my $chrStartEndG = $chrG .":". $startG ."-". $endG; #same as coords but without +/-
		my $chrStartEndA = $chrA .":". $startA ."-". $endA; 
		
		## check if passes filter ##
		if ($filter){
			my $found=0; 
			if($f_pairs){
				$found++ if (exists $f_pairs->{$coordsG ."\t". $coordsA} or exists $f_pairs->{$chrStartEndG ."\t". $chrStartEndA});
			}
			if($f_seqsG){
				$found++ if (exists $f_seqsG->{$coordsG} or exists $f_seqsG->{$chrStartEndG}); #with or without strand specifier
			}
			if($f_seqsA){
				$found++ if (exists $f_seqsA->{$coordsA} or exists $f_seqsA->{$chrStartEndA}); 
			}
			if($f_subfams){
				$found++ if exists $f_subfams->{$name}; #this is generic (it is possible that the hash will have only one subfam).	
			}
			
			#skip/don't skip this cluster by filter
			if (not $f_negate){ #if flag was 0 (remove arged pairs/seqs/subfams)
				next if $found; 
			}
			else{ #flag was 1 (retain only arged pairs/seqs/subfams)
				if($f_seqsG and $f_seqsA){ #need both seqs to exist to retain cluster
					next if $found < 2; #one or none of G and A seqs found
				}
				else{
					next if $found == 0; 
				}
			}
		}
		## print tabular cluster to file ##
		$retained++; #count retained cluster
		print TAB "$tuple\n";
		
		## Save cluster data in DSs ##
		$coords_G_all{$coordsG} = 0;
		$coords_A_all{$coordsA} = 0;
		foreach my $site (split(/,/,$locG)){
			$sites_G{$coordsG}{$site} = 0; 
		}
		foreach my $site (split(/,/,$locA)){
			$sites_A{$coordsA}{$site} = 0; 
		}
		
		$graphEntries{$coordsG."=$name"."\t".$coordsA."=$name"} = 0 unless exists $graphEntries{$coordsG."=$name"."\t".$coordsA."=$name"};
		($coordsToDefline{$coordsG}) = $g_name_l =~ /G name = (\S+)/;
		($coordsToDefline{$coordsA}) = $a_name_l =~ /A name = (\S+)/;
	}
}
close (TAB); 
close (CLUSTS); 
if ($filter and ($retained == 0)){ #exit if no clusters passed filter
	die "No clusters retained after filter for: filter: $filter and output dir: $dir\n"; 
}

### create track and graph files ###
open( my $seq_G_fh,    ">" . $dir . "/seq_G_".$suffix.".gff" );
open( my $seq_A_fh,    ">" . $dir . "/seq_A_".$suffix.".gff" );
open( my $sites_G_fh, ">" . $dir . "/sites_G_".$suffix.".gff" );
open( my $sites_A_fh, ">" . $dir . "/sites_A_".$suffix.".gff" );
my $siteListName_G = $dir . "/siteList_G_".$suffix.".txt"; 
open( my $siteList_G_fh, ">" . $siteListName_G );
my $siteListName_A = $dir . "/siteList_A_".$suffix.".txt"; 
open( my $siteList_A_fh, ">" . $siteListName_A );
open( my $graph_fh, ">" . $dir . "/graph_".$suffix.".txt");
open( my $graph2_fh, ">" . $dir . "/graph2_".$suffix.".txt");
open( my $outDeg_fh, ">" . $dir . "/outDeg_".$suffix.".txt");
open( my $inDeg_fh, ">" . $dir . "/inDeg_".$suffix.".txt");
#create header for GFF files. G sequences and sites: forest green; A sequences and sites: red.
my $desc = $organism . " " . $class . " " . $family . " " . $pval . " " . $th;
print $seq_G_fh "track name=\"Edited Sequences (G) " . $desc . "\" color=34,139,34 visibility=1\n"; 
print $seq_A_fh "track name=\"Edited Sequences (A) " . $desc . "\" color=255,0,0 visibility=1\n";
print $sites_G_fh "track name=\"Edited Sites (G) " . $desc . "\" color=34,139,34 visibility=1\n";
print $sites_A_fh "track name=\"Edited Sites (A) " . $desc . "\" color=255,0,0 visibility=1\n";

#print data to files
printCoordHashToGFF($seq_G_fh, \%coords_G_all); 
printCoordHashToGFF($seq_A_fh, \%coords_A_all); 
foreach my $entry (sort keys %graphEntries){
	print $graph_fh $entry ."\n";
	$entry =~ /([^=]+)=(\S+)\t([^=]+)=\S+/; 
	print $graph2_fh $organism ."\t". $class ."\t". $family ."\t". $2 ."\t". $1 ."\t". $3 ."\n"; 
}
printSitesToGFF($sites_G_fh, \%sites_G); 
printSitesToGFF($sites_A_fh, \%sites_A); 
printSiteList($siteList_G_fh, \%sites_G, \%coordsToDefline);
printSiteList($siteList_A_fh, \%sites_A, \%coordsToDefline); 
printInAndOutDeg($outDeg_fh, $inDeg_fh, \%graphEntries); 

### close files ### 
close($seq_G_fh); 
close($seq_A_fh); 
close($sites_G_fh); 
close($sites_A_fh); 
close($siteList_G_fh); 
close($siteList_A_fh); 
close($graph_fh); 
close($graph2_fh); 
close($outDeg_fh); 
close($inDeg_fh); 

### create fasta files for G and A sequences ###
my $seqFa_G_file = $dir . "/seqFasta_G_".$suffix.".fa";
my $seqFa_A_file = $dir . "/seqFasta_A_".$suffix.".fa"; 
system ("perl516 analysis_scripts/getFastaFromCoords.pl $siteListName_G $seqFa_G_file 0 0 0"); #***temp comment
system ("perl516 analysis_scripts/getFastaFromCoords.pl $siteListName_A $seqFa_A_file 0 0 0"); #***temp comment

### Find Context Preference Motifs of edited sites ###
system("perl516 analysis_scripts/findMotifs.pl $organism $class $family $pval $th $pmotif $control $subdir G $dir"); #***temp comment
system("perl516 analysis_scripts/findMotifs.pl $organism $class $family $pval $th $pmotif $control $subdir A $dir"); #***temp comment

### Advanced analyses - analysisSubs ###
getNucAndTrioFreq($dir, "G", $siteListName_G); 
getNucAndTrioFreq($dir, "A", $siteListName_A); 

# analysisSubs::getEditedPositionsInCons($siteListName_G); 
# analysisSubs::getEditedPositionsInCons($siteListName_A); 

### prints all keys in a hash to the arg file handle ### 
sub printCoordHashToGFF{
	(my $fh, my $hashRef) = @_; 
	my $count = 1; 
	foreach my $key (sort keys %$hashRef){
		# print $key ."\n"; #***
		(my $chr, my $start, my $end, my $strand) = $key =~ /^([^:]+):(\d+)-(\d+)([+-])/; 
		$start++; #increment to 1-base start for gff format
		print $fh $chr."\t"."DNAEfinder"."\t"."."."\t".$start."\t".$end."\t"."."."\t".$strand."\t"."."."\t".$count."\n";
		$count++; 
	}
}

#create two files - one with list of inDegs and one with outDegs
#fields: Name(subfamily)	Coords	InDeg/OutDeg
sub printInAndOutDeg{
	(my $out_fh, my $in_fh, my $graphEntryRef) = @_; 
	my $g_seq; my $a_seq; 
	my %out = (); my %in = (); 
	#count in degree of A sequences and out degree of G sequences
	foreach my $entry (keys %$graphEntryRef){
		($g_seq, $a_seq) = split(/\t/,$entry); 
		$out{$g_seq} = 0 unless exists $out{$g_seq};
		$out{$g_seq}++;
		$in{$a_seq} = 0 unless exists $in{$a_seq};
		$in{$a_seq}++;
	} 
	foreach my $g_seq (sort {$out{$b} <=> $out{$a}} keys %out){
		print $out_fh $g_seq."\t".$out{$g_seq}."\n"; 
	}
	foreach my $a_seq (sort {$in{$b} <=> $in{$a}} keys %in){
		print $in_fh $a_seq."\t".$in{$a_seq}."\n"; 
	}
}

sub printSitesToGFF{
	(my $fh, my $hashRef) = @_; 
	foreach my $coord (sort keys %$hashRef){
		(my $chr, my $start, my $end, my $strand) = $coord =~ /^([^:]+):(\d+)-(\d+)([+-])/; 
		my $site_pos;
		foreach my $site (sort {$a <=> $b} keys %{$hashRef->{$coord}}){
			if ($strand eq '+'){ #in positive strand - calculate from start ("start" coord is 0-base thus no need to decrease 1).
				$site_pos = $start + $site; 
			}
			else{ #in negative strand - calculate from end (decrease 1 because "end" coord is 1-base).
				$site_pos = $end - $site + 1; 
			}
			print $fh $chr."\t"."DNAEfinder"."\t"."."."\t".$site_pos."\t".$site_pos."\t"."."."\t".$strand."\t"."."."\t".$coord."\n";
		}
	}
}

sub printSiteList{
	(my $fh, my $hashRef, my $coordsToDefline) = @_;
	foreach my $coord (sort keys %$hashRef){
		print $fh $coordsToDefline->{$coord};
		foreach my $site (sort {$a <=> $b} keys %{$hashRef->{$coord}}){
			print $fh "\t" . $site;
		}
		print $fh "\n";
	}
}

sub createFilterFile{
	(my $dir, my $filter) = @_; 
	my $filter_file = $dir ."/". "filter.txt"; 
	open (FILTER, ">$filter_file") || die "open $filter_file\n"; 
	print FILTER "Filter=".$filter."\n"; 
	close(FILTER); 
}

#***Note: make sure that the return values are interpreted correctly for hashes that are empty
sub parseFilter{
	(my $filter) = @_;
	# print $filter ."\n"; 
	(my $type, my $negate, my $file) = split(/\|/,$filter); #(note: description was already stripped from filter)
	#print $type ."\t". $negate ."\t" . $file ."\n"; exit; 
	my $pairs=0; 
	my $seqsG=0; 
	my $seqsA=0; 
	my $subfams=0;
	my $file_lines; 
	if ($type ne "subfam"){ #get input from file
		$file_lines = getAll::lines($file, 1); #subfam doesn't input a file
		die "bad input file $file\n" unless $file_lines;
		
		#extract coords from full defline in format like: 81056=hg19=chr15:43315506-43315824+=SINE=Alu=AluSx1
		if ($file_lines->[0] =~ /=/){ 
			if ($type eq "pairs"){ #get coords from pairs
				for my $i (0 .. (scalar(@$file_lines) - 1)){
					(my $s, my $t) = $file_lines->[$i] =~ /=([^=\s]+:\d+-\d+[+-])=\S+\t\S+=([^=\s]+:\d+-\d+[+-])=/; 
					$file_lines->[$i] = $s ."\t". $t; 
				}
			}
			elsif($type =~ /seqs/){ #get sequences (for all options except for pairs, subfam and subfams)
				for my $i (0 .. (scalar(@$file_lines) - 1)){
					($file_lines->[$i]) = $file_lines->[$i] =~ /=([^=\s]+:\d+-\d+[+-])=/; 
				}
			}
		}
	}
	#save input from file to hash
	if ($type eq "pairs"){
		$pairs = {}; 
		foreach my $pair (@$file_lines){
			$pairs->{$pair}=0; 
		}
	}
	elsif($type eq "seqs"){
		$seqsG = {}; 
		$seqsA = {}; 
		foreach my $seq (@$file_lines){
			$seqsG->{$seq}=0; 
			$seqsA->{$seq}=0; 
		}
	}
	elsif($type eq "seqsG"){
		$seqsG = {}; 
		foreach my $seq (@$file_lines){
			$seqsG->{$seq}=0; 
		}
	}
	elsif($type eq "seqsA"){
		$seqsA = {}; 
		foreach my $seq (@$file_lines){
			$seqsA->{$seq}=0; 
		}
	}
	elsif($type eq "subfam"){
		$subfams = {}; 
		my $subfam = $file; #the file is actually the subfam
		$subfams->{$subfam}=0; 
	}
	elsif($type eq "subfams"){
		$subfams = {}; 
		foreach my $sf (@$file_lines){
			$subfams->{$sf}=0; 
		}
	}
	else{
		die "die: Bad filter type: $type - must be pairs, seqs or subfam\n"; 
	}
	return ($negate, $pairs, $seqsG, $seqsA, $subfams); 
}


### ADVANCED ANALYSES ###
sub getNucAndTrioFreq{ 
	(my $dir, my $GA, my $siteListFile) = @_; 
	my $suffix = $siteListFile;
	$suffix =~ s/.*siteList_//; 
	
	$dir .= "/"; 
	$GA = "G" unless $GA; #default "G"
	my $ga = lc $GA; 
	my $trim = 0; #if to trim polyA tail using trimest
	my $border = 0; #if to calc nuc freq only within border of first and last editing sites
	my $nucleotides = 0; #"cgt"; 
	
	# my $seqFile = $dir . "seqFasta_".$GA."_Human_SINE_Alu_1e-0_4.fa";
	my $seqFile = $dir . "seqFasta_". $suffix;
	$seqFile =~ s/.txt$/.fa/; 
	if ($trim){
		my $trimmedSeqFile = $seqFile; 
		$trimmedSeqFile =~ s/seqFasta/seqFastaTrimmed/; 
		my $mismatches = 1; 
		analysisSubs::trimPolyA($seqFile, $trimmedSeqFile, 0, $mismatches);
		$seqFile = $trimmedSeqFile; 
	}
	
	$siteListFile = $dir . "siteList_".$suffix; 
	my $freqFile = $dir . "logo_".$suffix; 
	my $range = 2; 
	analysisSubs::getNucFrequencyPerPosAllSeqs($seqFile, $siteListFile, $range, 1, 0, $freqFile); #(my $seqFile, my $siteListFile, my $range, my $normalize, my $reverse_editing, my $outfile)
	
	my $hash_ref; 
	if ($border){
		$hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $seqFile, $siteListFile, $nucleotides); #norm within borders!!!
	}
	else{
		$hash_ref = analysisSubs::normalizeFreqByFasta($freqFile, $seqFile, 0, $nucleotides);
	}
	 my $freqNormedFile = $siteListFile; 
	 $freqNormedFile =~ s/.*siteList/freq/; 
	 $freqNormedFile = $dir .($border ? "bordered_" : "").($nucleotides ? $nucleotides ."_" : ""). $freqNormedFile; 
	# my $freqNormedFile = $dir .($border ? "bordered_" : "").($nucleotides ? $nucleotides ."_" : ""). "freq_".$GA."_Human_SINE_Alu_1e-0_4.txt"; 
	analysisSubs::printFreqHash($hash_ref, $freqNormedFile);

	my $acgt = $ga; my $normalize = $seqFile; my $suppressPrint = 1; my $revcom = 0; my $range=1; my $reverse_editing=0;  
	# my $acgt = $ga; my $normalize = $seqFile; my $suppressPrint = 0; my $revcom = 0; my $range=1; my $reverse_editing=0; my $border=1; #norm within borders!!!
	my $triplets_ref = analysisSubs::getTripletsAllSeqs($seqFile, $acgt, $siteListFile, $range, $normalize, $reverse_editing, $suppressPrint, $revcom, $border);
	
}

