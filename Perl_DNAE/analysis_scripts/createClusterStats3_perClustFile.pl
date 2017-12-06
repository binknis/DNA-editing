#Function: parses a specific cluster file and creates clust_stats3 files for families and subfamilies (The file should contain all types of mismatches)
#	Input: 	1. clustFile - input cluster file
#			2. outPath - full path to out dir
#			3. suffix - description of file, to append to output files
#			4. f_or_sf - what output to print: 0 - both, 1 - fam, 2 - subfam
#			5. assemblyToOrgFile - enables conversion from assembly to org for output
#			6. constTaxaLabel - enables to input cluster file that doesn't org-class-fam-sf cols, but has a project definition col.
#	Output: 1. "cluster_stats3_fams_..." - Stat file for the class (regardless if any editing was found)
#			2. "cluster_stats3_subfams_..." - Stats for subfamilies that had at least some editing. 
#		Output fields: org	class	fam	subfam	pval th	edit_pairs	edit_pairs_c	FP_pairs	edited	edited_c	FP_edited	
#							----> edited_sites	edited_sites_c	FP_edited_sites	source	source_c	FP_source source_sites	source_sites_c FP_source_sites

use strict;
my $home = $ENV{"HOME"}; 
use lib  $home . "/Perl_DNAE"; 
use File::Path qw(mkpath);
(my $clustFile, my $outPath, my $suffix, my $f_or_sf, my $assemblyToOrgFile, my $constTaxaLabel) = @ARGV; 
my @nucPairs = ("GA", "CT", "GC", "GT", "CA", "TA"); 
my $famStatsFile = $outPath."/cluster_stats3_fams_".$suffix.".txt";
my $subfamStatsFile = $outPath."/cluster_stats3_subfams_".$suffix.".txt"; 
my $subfamStatsFileForR = $outPath."/cluster_stats3_subfams_".$suffix."_forR.txt"; 
mkpath $outPath; 
#Get Assembly to Organism map. 
#fields: clade org assembly

$assemblyToOrgFile = "/home/alu/binknis/Analyses/Utils/Map_Clades_Orgs_Assembly.txt" unless ($assemblyToOrgFile or $assemblyToOrgFile eq '0'); 
my %a2o = (); 
if(-e $assemblyToOrgFile){
	open (ATO, $assemblyToOrgFile) or die "$assemblyToOrgFile\n"; 
	while (my $l = <ATO>){
		chomp $l; 
		my @f = split(/\t/, $l); 
		$a2o{$f[2]} = $f[1]; 
	}
	close(ATO); 
}

my %fam_pairs = (); 
my %fam_seqs = (); 
my %fam_sites = (); 
my %sf_pairs = (); 
my %sf_seqs = (); 
my %sf_sites = (); 

#Cluster file fields: "[0] mismatch [1]assembly [2]class [3]family [4]subfam [5]coordsS [6]coordsT [7]num_this_mm [8]num_all_mms [9]total_prob [10]whereS [11]whereT
#						[12]num_clusts [13]mmSerials [14]clusters_span_woGaps [15]clusters_span [16]alignment_len   mmCount_str org      pval"
my $org; my $class; my $fam; my $subfam; my $coordsS; my $coordsT; my $whereS; my $whereT; 
my $taxa_f; my $taxa_sf; 
open(CLUSTS, $clustFile) or die "open $clustFile\n"; 
while (my $l = <CLUSTS>){
	chomp $l; 
	next if $l =~ /^mismatch/; #skip header
	
	my @f = split(/\t/, $l); 
	my $mm = $f[0]; 
	unless($constTaxaLabel){
		$org = $a2o{$f[1]};
		$class = $f[2]; 
		$fam = $f[3]; 
		$subfam = $f[4]; 
		$coordsS = $f[5]; 
		$coordsT = $f[6]; 
		$whereS = $f[10]; 
		$whereT = $f[11]; 
		
		$taxa_f = join("\t", $org, $class, $fam); 
		$taxa_sf = join("\t", $org, $class, $fam, $subfam); 
	} else {
		$taxa_f = $taxa_sf = $constTaxaLabel; 
		$coordsS = $f[2]; 
		$coordsT = $f[3]; 
		$whereS = $f[7]; 
		$whereT = $f[8]; 
		
	}

	

	##save pair 
	$fam_pairs{$taxa_f}{$mm}{$coordsS ."\t". $coordsT}=0;
	$sf_pairs{$taxa_sf}{$mm}{$coordsS ."\t". $coordsT}=0; 
	
	##save elements
	$fam_seqs{$taxa_f}{$mm}{"s"}{$coordsS}=0;
	$fam_seqs{$taxa_f}{$mm}{"t"}{$coordsT}=0;
	$sf_seqs{$taxa_sf}{$mm}{"s"}{$coordsS}=0;
	$sf_seqs{$taxa_sf}{$mm}{"t"}{$coordsT}=0;
	
	##save sites
	#source sites
	(my $chr, my $start, my $end, my $strand) = $coordsS =~ /^([^:]+):(\d+)-(\d+)([+-])$/; 
	my @s_sites = split(',', $whereS); 
	foreach my $pos (@s_sites){
		my $siteCoord = $chr . ":" . ($strand eq '+' ? ($start + $pos) : ($end - $pos + 1) ); #get correct position based on + or - 
		$fam_sites{$taxa_f}{$mm}{"s"}{$siteCoord} = 0; 
		$sf_sites{$taxa_sf}{$mm}{"s"}{$siteCoord} = 0; 
	}
	#target sites
	(my $chr, my $start, my $end, my $strand) = $coordsT =~ /^([^:]+):(\d+)-(\d+)([+-])$/; 
	my @t_sites = split(',', $whereT);
	foreach my $pos (@t_sites){
		my $siteCoord = $chr . ":" . ($strand eq '+' ? ($start + $pos) : ($end - $pos + 1) ); #get correct position based on + or - 
		$fam_sites{$taxa_f}{$mm}{"t"}{$siteCoord} = 0;
		$sf_sites{$taxa_sf}{$mm}{"t"}{$siteCoord} = 0;
	}
}
close(CLUSTS); 

## print output for all families
if ($f_or_sf==0 or $f_or_sf==1){
	open (my $f_fh, ">".$famStatsFile) or die "open $famStatsFile\n"; 
	
	#print header
	my $header; 
	unless($constTaxaLabel){
		$header = join("\t", "org", "class", "fam"); 
	} else {
		$header = join("\t", "Project"); 
	}
	foreach my $mm (@nucPairs){
		foreach my $suff ("pairs", "s_elements", "t_elements", "s_sites", "t_sites"){
			$header .= "\t" . $mm . "_" . $suff; 
		}
	}
	print $f_fh $header . "\n"; 
	
	#print results
	foreach my $taxa_f (sort keys %fam_pairs){
		print $f_fh $taxa_f; 
		foreach my $mm (@nucPairs){
			my $num_pairs = keys %{$fam_pairs{$taxa_f}{$mm}};
			my $num_sseqs = keys %{$fam_seqs{$taxa_f}{$mm}{"s"}};
			my $num_tseqs = keys %{$fam_seqs{$taxa_f}{$mm}{"t"}};
			my $num_ssites = keys %{$fam_sites{$taxa_f}{$mm}{"s"}};
			my $num_tsites = keys %{$fam_sites{$taxa_f}{$mm}{"t"}};
			my @nums = ($num_pairs, $num_sseqs, $num_tseqs, $num_ssites, $num_tsites); 
			for my $i(0 .. $#nums){ #convert empty to zeros
				$nums[$i] = ($nums[$i] ? $nums[$i] : 0); 
			}
			print $f_fh "\t" , join("\t", @nums); 
		}
		print $f_fh "\n"; 
	}
	close($f_fh);
}

## print output for all subfamilies
if (($f_or_sf==0 or $f_or_sf==2) and not $constTaxaLabel){
	open (my $sf_fh, ">".$subfamStatsFile) or die "open $subfamStatsFile\n"; 
	open (my $sf_fh_R, ">".$subfamStatsFileForR) or die "open $subfamStatsFileForR\n"; 
	#create and print header
	my $header = join("\t", "org", "class", "family", "subfam"); 
	my @cols = ("pairs", "s_elements", "t_elements", "s_sites", "t_sites"); 
	#print short header - for R 
	print $sf_fh_R join("\t", "mismatch", $header, @cols) , "\n"; 
	#continue creating header for regular file 
	foreach my $mm (@nucPairs){
		foreach my $suff (@cols){
			$header .= "\t" . $mm . "_" . $suff; 
		}
	}
	print $sf_fh $header . "\n"; 

	#print results
	foreach my $taxa_sf (sort keys %sf_pairs){
		print $sf_fh $taxa_sf; 
		foreach my $mm (@nucPairs){
			my $num_pairs = keys %{$sf_pairs{$taxa_sf}{$mm}};
			my $num_sseqs = keys %{$sf_seqs{$taxa_sf}{$mm}{"s"}};
			my $num_tseqs = keys %{$sf_seqs{$taxa_sf}{$mm}{"t"}};
			my $num_ssites = keys %{$sf_sites{$taxa_sf}{$mm}{"s"}};
			my $num_tsites = keys %{$sf_sites{$taxa_sf}{$mm}{"t"}};
			my @nums = ($num_pairs, $num_sseqs, $num_tseqs, $num_ssites, $num_tsites); 
			for my $i(0 .. $#nums){ #convert empty to zeros
				$nums[$i] = ($nums[$i] ? $nums[$i] : 0); 
			}
			print $sf_fh "\t" , join("\t", @nums); 
			print $sf_fh_R join("\t", $mm, $taxa_sf, @nums) , "\n"; 
		}
		print $sf_fh "\n"; 
	}
	close($sf_fh);
	close($sf_fh_R); 
}
