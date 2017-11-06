#Function: Parses a specific cluster file and creates several UCSC-track and analysis-output files
#Input: Specification of a cluster file (see args)
#		subdir = in result file
#Output: creates - 	1. Tabular cluster file; 
#					2. graph file (0-base starts); 
#					3. graph2 - fields: org class fam name assembly source_coords target_coords (Note: assumes name of source and target are the same)
#					3. track files (gff files are 1-base start); 
#					4. 
#					5. 
#exec-example:  perl analysis_scripts/createTrackFiles.pl Human LINE L1 1e-6 10  0 Clusters_80_percent_length_G &
use strict;
( my $organism, my $class, my $family, my $pval, my $th, my $control, my $subdir) = @ARGV;
$pval = "1e-" . $pval if $pval =~ /^\d+$/;
$subdir = 0 if $subdir eq ''; 
$control = 0 if $control eq '';
my $pmotif = '1e-3'; 
#### build file name and track dir ####
#file name 

my $cluster_file = "../Data/" . $organism . "/" . $class . "/results"; 
$cluster_file .= "/$subdir" unless ($subdir eq 0 or $subdir eq ""); 
$cluster_file .= "/clusters_" . $organism . "_" . $class . "_" . $family . "_" . $pval . "_" . $th;
#track dir
my $dir = "../Data/" . $organism . "/" . $class . "/results"; 
$dir .= "/$subdir" unless ($subdir eq 0 or $subdir eq ""); 
$dir .= "/Tracks"; 
mkdir $dir; #make "Tracks" dir; later will create specific 'tracks' dir
$dir .= "/tracks_" . $organism . "_" . $class . "_" . $family . "_" . $pval . "_" . $th;
#file and dir trailers
if ($control)
{
	$cluster_file .= "_control";
	$dir     .= "_control";
}

(my $suffix) = $cluster_file =~ /clusters_(\S+)$/;
my $tabular_file = $dir ."/clusters_". $suffix .".tab"; 
$cluster_file .= ".txt";
#don't parse if cluster's file is empty (avoids creating empty track files). 
exit if (-z $cluster_file); 
mkdir $dir;

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
open (CLUSTS ,"<" . $cluster_file) || die ("couldn't open $cluster_file"); 
open (TAB ,">" . $tabular_file) || die ("couldn't open $tabular_file"); 
my $tuple;
while(my $line  = <CLUSTS>)
{
	if($line =~ /^Found/)
	{
		### READ CLUSTER ### 
		my $g_name_l  = <CLUSTS>;
		chomp $g_name_l; 
		my $a_name_l = <CLUSTS>;
		chomp $a_name_l; 
		my $prob_clustNum_l = <CLUSTS>;
		chomp $prob_clustNum_l; 
		my $mm_serial_num_l = <CLUSTS>;
		chomp $mm_serial_num_l; 
		my $locG_l = <CLUSTS>;
		chomp $locG_l; 
		my $locA_l = <CLUSTS>;
		chomp $locA_l; 
		my $len_mmNum_l = <CLUSTS>;
		chomp $len_mmNum_l; 
	
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
		# print $data "$tuple\n";
		#print tabular cluster to file
		print TAB "$tuple\n";
		
		##### Save cluster data in DSs #####
		$coordsG = $chrG .":". $startG ."-". $endG .$strandG; 
		$coordsA = $chrA .":". $startA ."-". $endA .$strandA; 
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
system ("perl516 analysis_scripts/getFastaFromCoords.pl $siteListName_G $seqFa_G_file 0 0 0"); 
system ("perl516 analysis_scripts/getFastaFromCoords.pl $siteListName_A $seqFa_A_file 0 0 0"); 

### Find Context Preference Motifs of edited sites ###
system("perl516 analysis_scripts/findMotifs.pl $organism $class $family $pval $th $pmotif $control $subdir G"); 
system("perl516 analysis_scripts/findMotifs.pl $organism $class $family $pval $th $pmotif $control $subdir A"); 


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
