#Function: Parses a specific cluster file and creates several UCSC-track and analysis-output files
#			CHANGE FROM createTrackFiles.pl: Works on tabular output cluster files
#Input: Specification of a cluster file (see args)
#	filter = select only part of sequences in cluster file: Format := [pairs/seqs/seqsG/seqsA/subfam/subfams]|[1/0]|[bed-file/gff-file/pair-file/subfam/subfam-file]|description. example: "seqsA|1|file.bed|this_describes_filter" (means: a sequence file for A sequences | retain only these | BED format input | "this_describes_filter" will be appended to output dir)
#     Filter explanation and pairing: 
	#		1. pairs - 3rd param will be two columns in "chrN:start-end[+-]" format 
	#		2. seqs - filter for both G and A seqs (BED or GFF file)
	#		3. seqsG - filter for G seqs 
	#		4. seqsA - filter for A seqs 
	#		5. subfam - one subfam
	#		6. subfams - a file with a list of subfams		
	#		7. clusters - list of clusters to create tracks from (reads directly from the clusters file specified)
	#					clusters-file spcifications: must contain cluster for this org and class only, but can contain for multiple families (only clusters for arged family will be analyzed)
	#   1/0: retain (1) or remove (0) the specified sequences
	# Notes: 
	#	1. Input format of sequences in filter-file must be full(with multiple "=" signs) or simply coords (with or without strand)
	#	2. For "clusters" filter, the 'retain' flag is disregarded and set to 1. 
#		Filters already completed: Maybe all... need to verify ***
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
use lib $ENV{HOME} ."/Perl_DNAE";
use getAll; 
my $perlDir = $ENV{HOME} ."/Perl_DNAE"; 
use lib $ENV{HOME} ."/Perl_DNAE/analysis_scripts"; 
use analysisSubs; 
use File::Path qw(mkpath);
my $home = $ENV{HOME}; 
# use Data::Dumper;

# my $INTERVAL_DIR = "/private3/Dropbox/users/binknis_data/Avian/rawdata/intervals"; #interval dir (currently for Avian genome project)
my $INTERVAL_DIR = "/home/alu/binknis/binknis_data/Raw_Data/Vertebrate"; #interval dir (currently for Avian genome project)
my $FETCH_SUBFAMS_FROM_INTERVAL_FILE = 0; #CONST to get or not to get subfams from the original Interval file  (1 was used for avian genomes)
my $SUBFAM_SELECTION_NUC; #if there are different subfams - use G or A's subfam when single annotation is needed
my $DONT_GET_BEST_SOURCES = 1; #skip last time-consuming step that BLASTs to find best match (good to disable for large intermediate datasets)
( my $org, my $class, my $family, my $pval, my $th, my $mm, my $filter) = @ARGV;
$pval = "1e-" . $pval if $pval =~ /^\d+$/;
$mm = "GA" if ($mm eq '' or $mm eq '0');
$mm = uc $mm; 
(my $mmS, my $mmT) = split('', $mm); #save source and target mismatches
unless ($SUBFAM_SELECTION_NUC){ 
	$SUBFAM_SELECTION_NUC = $mmS; 
}
my $pmotif = '1e-3';
my $pmotif2 = '1e-2';
my $filter_subdir = ($filter ? "/Filtered" : ""); 
my $f_type;  my $retain; my $f_pairs=""; my $f_seqsG=""; my $f_seqsA=""; my $f_subfams=""; my $cluster_file_filter; #filter variables used later
my %subfamPerElem;
#### build file name and track dir ####
#cluster file
my $suffix = $org . "_" . $class . "_" . $family . "_" . $pval . "_" . $th; 
my $resDir = $home . "/Data/" . $org . "/" . $class . "/results"; 
my $cluster_file = $resDir ."/". $mm . "/clusters_" . $suffix .".tab";
exit if (not -e $cluster_file or -z $cluster_file); #don't parse if cluster's file is empty (avoids creating empty track files). 
#track dir ("dir")
my $Tracks_dir = $resDir . "/Tracks"; 
my $dir = $Tracks_dir . $filter_subdir; 
$dir .= "/tracks_" . $suffix ."/". $mm;
#dir trailers
if ($filter){
	(my $f_suffix) = $filter =~ /([^\|]+)$/;
	$dir .= "/" . $f_suffix;
	$filter =~ s/\|[^\|]+$//; 
}
#create track dirs
mkdir $Tracks_dir; #create "Tracks" dir
mkpath $dir; #create specific 'tracks' dir
my $tabular_file = $dir ."/clusters_". $suffix .".tab"; #cluster_file_tabular_in_track_dir

#create documentation file for filter and parse the filter (if file is specified it will be parsed)
if ($filter){
	createFilterFile($dir, $filter); #Create filter-description file
	($retain, $f_pairs, $f_seqsG, $f_seqsA, $f_subfams, $cluster_file_filter) = parseFilter($filter);
	if ($cluster_file_filter){
		$cluster_file = $cluster_file_filter;
	}
}

## DSs for graph and tracks ## 
my %coords_G_all  = ();
my %coords_A_all  = (); 
my %minProb_G = (); 
my %minProb_A = (); 
my %sites_G = (); 
my %sites_A = (); 
my %graphEntries = ();
my %coordsToDefline = (); 

my $coordsToSubfam = (); 
if($FETCH_SUBFAMS_FROM_INTERVAL_FILE){
	$coordsToSubfam = getSubfamsFromIntervalFile($INTERVAL_DIR, $org, $class, $family); 
}

#### Parse cluster file and create tabular cluster file ###
my $retained = 0; #for filter count (to exit of 0 are retained after filter)
open (CLUSTS ,"<" . $cluster_file) || die ("couldn't open $cluster_file"); 
if ($filter or $FETCH_SUBFAMS_FROM_INTERVAL_FILE){ #For filter - create a new file with the clusters retained for creating files in trackDir
	open (TAB ,">" . $tabular_file) || die ("couldn't open $tabular_file"); 
}
else{ #No filter - there is no need to create a new cluster file, a symbolic link will be created to original cluster file in the trackDir
	unlink($tabular_file); #in case it was already created and needs to override existing (needed after correcting bug).
	symlink($cluster_file, $tabular_file);
}
while(my $line  = <CLUSTS>)
{
	chomp $line; 
	next unless $line =~ /^[ACTG][ACTG]\t/; #skip header line
	# print $line; #***
	(my $mmType, my $assembly, my $class, my $fam, my $subfam, 
	my $coordsG, my $coordsA, my $direct_mms, my $num_all_mms, my $prob, 
	my $locG, my $locA, my $num_clusts, my $mmSerials, 
	my $clusters_span_woGaps, my $clusters_span, my $align_length, my $mmCount_str) = split ('\t', $line); 
	next unless ($family eq $fam); #This is needed only whein "clusters" filter is on and a file with multiple families is specified
	
	my $subfamG; my $subfamA; 
	if($FETCH_SUBFAMS_FROM_INTERVAL_FILE){ #different subfams may occur for G and A seqs
		$subfamG = $coordsToSubfam->{$coordsG}; 
		$subfamPerElem{$mmS} = $subfamG; 
		$subfamA = $coordsToSubfam->{$coordsA}; 
		$subfamPerElem{$mmT} = $subfamA; 
		$subfam = $subfamPerElem{$SUBFAM_SELECTION_NUC};
	}
	else{ #same subfam for both G and A sequences
		$subfamG = $subfam; 
		$subfamA = $subfam;
	}
	
	(my $chrStartEndG) = $coordsG =~ /(\S+:\d+-\d+)[+-]?/; #same as coords but without +/-
	(my $chrStartEndA) = $coordsA =~ /(\S+:\d+-\d+)[+-]?/; 
	
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
			$found++ if exists $f_subfams->{$subfam}; #this is generic (it is possible that the hash will have only one subfam).	
		}
		if($cluster_file_filter){ #retain all for clusters filter
			$found=1; 
		}
		
		#skip/don't skip this cluster by filter
		if (not $retain){ #if flag was 0 (remove arged pairs/seqs/subfams)
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
		
		## passed filter - print tabular cluster to file ## 
		my $tuple; 
		if($cluster_file_filter or $FETCH_SUBFAMS_FROM_INTERVAL_FILE){
			$tuple = join("\t", $mmType, $assembly, $class, $fam, $subfam, $coordsG, $coordsA, 
								$direct_mms, $num_all_mms, $prob, $locG, $locA, $num_clusts, 
								$mmSerials, $clusters_span_woGaps, $clusters_span, $align_length, $mmCount_str);  #This is same order as parsing from original clusters line, but only first columns are saved to avoid printing analysis columns from processed clusters file
		}
		else{
			$tuple = $line; 
		}
		print TAB "$tuple\n";
	}
	elsif($FETCH_SUBFAMS_FROM_INTERVAL_FILE){
		my $tuple = join("\t", $mmType, $assembly, $class, $fam, $subfam, $coordsG, $coordsA, 
								$direct_mms, $num_all_mms, $prob, $locG, $locA, $num_clusts, 
								$mmSerials, $clusters_span_woGaps, $clusters_span, $align_length, $mmCount_str);  #This is same order as parsing from original clusters line, but only first columns are saved to avoid printing analysis columns from processed clusters file
		print TAB "$tuple\n";
	}
	$retained++; #count retained cluster
	
	## Save cluster data in DSs ##
	$coords_G_all{$coordsG} = 0;
	$coords_A_all{$coordsA} = 0;
	foreach my $site (split(/,/,$locG)){
		$sites_G{$coordsG}{$site} = 0; 
	}
	foreach my $site (split(/,/,$locA)){
		$sites_A{$coordsA}{$site} = 0; 
	}
	
	$graphEntries{$coordsG."=$subfamG"."\t".$coordsA."=$subfamA"} = 0;
	($coordsToDefline{$coordsG}) = $assembly ."=". $coordsG ."=". $class ."=". $fam ."=". $subfamG;
	($coordsToDefline{$coordsA}) = $assembly ."=". $coordsA ."=". $class ."=". $fam ."=". $subfamA; 	
}
close (CLUSTS); 
close (TAB) if $filter; 

if ($filter and ($retained == 0)){ #exit if no clusters passed filter
	die "No clusters retained after filter for: filter: $filter and output dir: $dir\n"; 
}

### create track and graph files ###
open( my $seq_G_fh,    ">" . $dir . "/seq_".$mmS."_".$suffix.".gff" );
open( my $seq_A_fh,    ">" . $dir . "/seq_".$mmT."_".$suffix.".gff" );
open( my $sites_G_fh, ">" . $dir . "/sites_".$mmS."_".$suffix.".gff" );
open( my $sites_A_fh, ">" . $dir . "/sites_".$mmT."_".$suffix.".gff" );
my $siteListName_G = $dir . "/siteList_".$mmS."_".$suffix.".txt"; 
open( my $siteList_G_fh, ">" . $siteListName_G );
my $siteListName_A = $dir . "/siteList_".$mmT."_".$suffix.".txt"; 
open( my $siteList_A_fh, ">" . $siteListName_A );
my $graphFile = $dir . "/graph_".$suffix.".txt"; 
open( my $graph_fh, ">" . $graphFile);
open( my $graph2_fh, ">" . $dir . "/graph2_".$suffix.".txt");
open( my $outDeg_fh, ">" . $dir . "/outDeg_".$suffix.".txt");
open( my $inDeg_fh, ">" . $dir . "/inDeg_".$suffix.".txt");
#create header for GFF files. G sequences and sites: forest green; A sequences and sites: red.
my $desc = $org . " " . $class . " " . $family . " " . $pval . " " . $th;
print $seq_G_fh "track name=\"Edited Sequences (".$mmS.") " . $desc . "\" color=34,139,34 visibility=1\n"; 
print $seq_A_fh "track name=\"Edited Sequences (".$mmT.") " . $desc . "\" color=255,0,0 visibility=1\n";
print $sites_G_fh "track name=\"Edited Sites (".$mmS.") " . $desc . "\" color=34,139,34 visibility=1\n";
print $sites_A_fh "track name=\"Edited Sites (".$mmT.") " . $desc . "\" color=255,0,0 visibility=1\n";

#print data to files
printCoordHashToGFF($seq_G_fh, \%coords_G_all); 
printCoordHashToGFF($seq_A_fh, \%coords_A_all); 
foreach my $entry (sort keys %graphEntries){
	print $graph_fh $entry ."\n";
	$entry =~ /([^=]+)=(\S+)\t([^=]+)=\S+/; #get coordsG, subfam, coordsA
	print $graph2_fh $org ."\t". $class ."\t". $family ."\t". $2 ."\t". $1 ."\t". $3 ."\n"; 
}
printSitesToGFF($sites_G_fh, \%sites_G); 
printSitesToGFF($sites_A_fh, \%sites_A); 
# print Dumper(\%coordsToDefline); exit; #***
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
my $seqFa_G_file = $dir . "/seqFasta_".$mmS."_".$suffix.".fa";
my $seqFa_A_file = $dir . "/seqFasta_".$mmT."_".$suffix.".fa"; 
unless($FETCH_SUBFAMS_FROM_INTERVAL_FILE){
	# system ("perl516 $perlDir/analysis_scripts/getFastaFromCoords.pl $siteListName_G $seqFa_G_file 0 0 0 0"); 
	# system ("perl516 $perlDir/analysis_scripts/getFastaFromCoords.pl $siteListName_A $seqFa_A_file 0 0 0 0"); 
	analysisSubs::getFastaFromCoords($siteListName_G, $seqFa_G_file, 0, 0, 0, 0);
	analysisSubs::getFastaFromCoords($siteListName_A, $seqFa_A_file, 0, 0, 0, 0);
}
else{
	analysisSubs::getFastaFromCoords($siteListName_G, $seqFa_G_file, 0, 0, 0, $coordsToSubfam);
	analysisSubs::getFastaFromCoords($siteListName_A, $seqFa_A_file, 0, 0, 0, $coordsToSubfam);
}

### Find Context Preference Motifs of edited sites ###
unless($FETCH_SUBFAMS_FROM_INTERVAL_FILE){ #The findMotifs.pl script wasn't adapted for this option. 
	system("perl516 $perlDir/analysis_scripts/findMotifs.pl $org $class $family $pval $th $pmotif $mm $mmS $dir");
	system("perl516 $perlDir/analysis_scripts/findMotifs.pl $org $class $family $pval $th $pmotif $mm $mmT $dir"); 

	if($pmotif2 > 0){
		system("perl516 $perlDir/analysis_scripts/findMotifs.pl $org $class $family $pval $th $pmotif2 $mm $mmS $dir");
		system("perl516 $perlDir/analysis_scripts/findMotifs.pl $org $class $family $pval $th $pmotif2 $mm $mmT $dir");
	}
}

### Advanced analyses - analysisSubs ###
analysisSubs::getEditedPositionsInCons($siteListName_G); 
analysisSubs::getEditedPositionsInCons($siteListName_A); 

getNucStats($dir, $mmS, $siteListName_G); 
getNucStats($dir, $mmT, $siteListName_A); 


#Get best sources for each target element
my $STreversed = 0; 
my $transform = 0; 
my $trimmed = 0; 
analysisSubs::getBestSources($graphFile, $seqFa_G_file, $seqFa_A_file, $STreversed, $transform, $trimmed) unless $DONT_GET_BEST_SOURCES; 

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
#***Note: make sure it works after change to all types of mismatches + tabular format
sub parseFilter{
	(my $filter) = @_;
	# print $filter ."\n"; 
	(my $type, my $retain, my $file) = split(/\|/,$filter); #(note: description was already stripped from filter)
	# print $type ."\t". $retain ."\t" . $file ."\n"; exit; #***
	my $pairs=0; 
	my $seqsG=0; 
	my $seqsA=0; 
	my $subfams=0;
	my $ret_cluster_file=0; 
	my $file_lines; 
	if ($type ne "subfam"){ #get input from file
		$file_lines = getAll::lines($file, 1); #subfam doesn't input a file
		die "bad input file $file\n" unless $file_lines;
		
		#extract coords from full defline in format like: 81056=hg19=chr15:43315506-43315824+=SINE=Alu=AluSx1
		if ($file_lines->[0] =~ /=/){ #*** check if this works after change to tabular format
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
	elsif($type eq "clusters"){
		$ret_cluster_file = $file; 
		$retain = 1; 
	}
	else{
		die "die: Bad filter type: $type - must be pairs, seqs, subfamG, subfamA or clusters\n"; 
	}
	return ($retain, $pairs, $seqsG, $seqsA, $subfams, $ret_cluster_file); 
}

### ADVANCED ANALYSES ###
sub getNucStats { 
	(my $dir, my $GA, my $siteListFile) = @_; 
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
	
	#create nucListFreq file
	my $nucListFile = $dir . "nucList_".$suffix.".txt"; 
	analysisSubs::nucListToFreq($nucListFile);
	
	#create nucList + nucListFreq per pairs
	(my $suffix2) = $suffix =~ /[ACTG]_(\S+)/;
	my $clustersFile = $dir . "clusters_".$suffix2.".tab"; 
	my $SorT = ($GA =~ /[AT]/ ? 1 : 0); 
	analysisSubs::nucFreqPerPairs($clustersFile, $SorT, $GA);
	
	#print nuc composition for all seqs and per seq
	my $ncFile = $seqFile; 
	$ncFile =~ s/seqFasta/nucComp/; 
	$ncFile =~ s/\.fa$/.txt/; 
	my $ncPerSeqFile = $seqFile; 
	$ncPerSeqFile =~ s/seqFasta/nucCompPerSeq/; 
	$ncPerSeqFile =~ s/\.fa$/.txt/;
	#Note: Fields for printAllSeqNucComposition (my $seqFile, my $retfmt, my $siteListFile, my $outFile, my $outFilePerSeq)
	analysisSubs::printAllSeqNucComposition($seqFile,1,0,$ncFile, $ncPerSeqFile); #no bordering 
	$ncFile =~ s/nucComp/nucComp_bordered/; 
	$ncPerSeqFile =~ s/nucCompPerSeq/nucCompPerSeq_bordered/; 
	analysisSubs::printAllSeqNucComposition($seqFile, 1, $siteListFile, $ncFile, $ncPerSeqFile); #composition between terminal editing sites
	
}

#***Add option to read from zipped file
sub getSubfamsFromIntervalFile{
	(my $invlDir, my $org, my $class, my $fam) = @_; 
	my %c2sf = (); 
	my $invlFile = $invlDir ."/". $org . ".interval"; 
	open(my $invl_fh, $invlFile) or die "open $invlFile\n"; 
	while(my $l = <$invl_fh>){
		chomp $l; 
		next if $l =~ /^#/; 
		my @fs = split(/\t/, $l); #chr, start, end, strand, name, class, fam 
		next unless ($fs[5] eq $class and $fs[6] eq $fam); #insert only this coords of this class and fam to hash
		my $coords = $fs[0] .":". $fs[1] ."-". $fs[2] . $fs[3]; 
		$c2sf{$coords} = $fs[4]; 
	}
	close($invl_fh); 
	return \%c2sf; 
}