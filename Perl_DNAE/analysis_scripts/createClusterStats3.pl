#Function: parses all cluster files of a class (of one org) and creates clust_stats3 files for families and subfamilies
#	Output: 1. "cluster_stats3_fams_..." - Stat file for the class (regardless if any editing was found)
#			2. "cluster_stats3_subfams_..." - Stats for subfamilies that had at least some editing. 
#		Output fields: Org	class	fam	subfam	pval th	edit_pairs	edit_pairs_c	FP_pairs	edited	edited_c	FP_edited	
#							----> edited_sites	edited_sites_c	FP_edited_sites	source	source_c	FP_source source_sites	source_sites_c FP_source_sites
#	Notes: 1. Removes all stat files currently existing (Family/Subfamily stats). 
#started changing. Find a way to save sites{G} and sites{A}
#			2. Names must have 0-base starts.  

use strict;
my $home = $ENV{"HOME"}; 
use lib  $home . "/Perl_DNAE"; 
use lib "/home/alu/binknis/Perl_DNAE"; 
use getAll;  
(my $org , my $class, my $paramLimits, my $subdir, my $dataDir) = @ARGV; 

$dataDir = $home ."/". "Data" unless $dataDir; 

my @nucPairs = ("GA", "CT", "GC", "GT", "CA", "TA"); 
my $resDir = join('/', $dataDir, $org, $class, "results"); 
my $famStatsFile = $resDir."/cluster_stats3_fams_".$org."_".$class.".txt"; 
my $subfamStatsFile = $resDir."/cluster_stats3_subfams_".$org."_".$class.".txt"; 
#remove old cluster_stats files
system("rm -f $famStatsFile");
system("rm -f $subfamStatsFile");

#get cluster files (full path)
my $pf = 1; #full path
my $clustFiles = getAll::tabClustersClass($org, $class, $pf, $subdir); 
die "no cluster_files for $org $class\n" if $clustFiles == 0; 
my $allParams = getParams($clustFiles, $paramLimits); 
my $fams = getAll::families($org, $class); 
my %fam_stats = ();
my %sf_stats = (); 

foreach my $fam (@$fams){ 
	my $subfams = getAll::names($org, $class, $fam);
	foreach my $np (@nucPairs){
		foreach my $cf (@{$clustFiles->{$np}}){
			(my $o, my $c, my $f, my $pval, my $th) = $cf =~ /clusters_([^_]+)_([^_]+)_(\S+)_(1e-\d+)_(\d+).tab$/;
			my $params = $pval ."\t". $th;
			# print "$pval\t$th\n"; 	
			# print Dumper($allParams); 
			next unless ($f eq $fam and exists $allParams->{$params});
			
			$sf_stats{$fam}{$params}{$np} = getSubfamStats($cf); 
			
			#init %fam_stats HoHoHoH
			$fam_stats{$fam}{$params}{$np}{"pairs"} = 0; 
			foreach my $s_or_t ("s", "t"){
				foreach my $s_or_s("seqs", "sites"){
					$fam_stats{$fam}{$params}{$np}{$s_or_t}{$s_or_s} = 0;
				}
			}
			#Sum stats for family
			foreach my $sf (@$subfams){
				# Sum data for family
				$fam_stats{$fam}{$params}{$np}{"pairs"} += $sf_stats{$fam}{$params}{$np}{$sf}{"pairs"};
				foreach my $s_or_t ("s", "t"){
					foreach my $s_or_s("seqs", "sites"){
						$fam_stats{$fam}{$params}{$np}{$s_or_t}{$s_or_s} += $sf_stats{$fam}{$params}{$np}{$sf}{$s_or_t}{$s_or_s};
					}
				}
			}
		}
	}
}
### print output for all families and subfamilies
open (my $f_fh, ">".$famStatsFile) or die "open $famStatsFile\n"; 
open (my $sf_fh, ">".$subfamStatsFile) or die "open $subfamStatsFile\n"; 
foreach my $fam (sort @$fams){
	printStats($f_fh, $org, $class, $fam, $allParams, \%{$fam_stats{$fam}}); #Print family cluster stats to family stat file
	my $subfams = getAll::names($org, $class, $fam);
	foreach my $sf (sort @$subfams){
		printStats($sf_fh, $org, $class, $fam, $allParams, \%{$sf_stats{$fam}}, $sf) ; #Print subfamily cluster stats to subfamily stat file
	}
}
close($f_fh);
close($sf_fh);

#print output (for family or subfamily; depends on $sf param)
sub printStats{
	(my $stats_fh, my $org, my $class, my $fam, my $allParams, my $stats, my $sf) = @_;
	foreach my $params (sort keys %$allParams){
		#build output tuple
		my $tuple = $org."\t".$class."\t".$fam."\t". ($sf ? $sf."\t" : '') .$params;
		foreach my $np (@nucPairs){
			#get correct stat_ref
			my $sr; 
			if ($sf){ #sf_stats
				$sr = $stats->{$params}{$np}{$sf}; 
			}
			else{ #fam_stats
				$sr = $stats->{$params}{$np}; 
			}
			#continue tuple build
			my $val = $sr->{"pairs"};  #use $val here and below to place 0 instead of blank output
			$tuple .= "\t" . ($val ? $val : '0'); 
			foreach my $s_or_t ("s", "t"){
				foreach my $s_or_s("seqs", "sites"){ 
					$val = $sr->{$s_or_t}{$s_or_s}; 
					$tuple .= "\t" . ($val ? $val : '0'); 
				}
			}
		}
		print $stats_fh $tuple . "\n";
	}
}

#Input: a clusters filename 
#Return: Hash sf_stats holds (for each subfam): a. pairs; b. source seqs; c. source sites; d. target seqs ; e. target sites. 
#Note: The hashes returned contain only the subfams that had clusters 
sub getSubfamStats{
	(my $cf_name) = @_; 
	my %sf_stats = (); #"pairs", {"s"}({seqs}, {sites}), {"t"}({seqs}, {sites}),
	my $subfam, my $seqs; my $sites;
	#return unless (-e $cf_fh); #return if file doesn't exist (for these params)
	open (my $cf_fh, "<". $cf_name) || die "didn't open $cf_name\n"; 
	do{
		($subfam, $seqs, $sites) = readTabCluster($cf_fh);
		if ($subfam){ #cluster fetched
			#add pair
			$sf_stats{$subfam}{"pairs"}{$seqs->{"s"} ."\t". $seqs->{"t"}}=0; 
			#Add seq and sites for source and target
			foreach my $s_or_t("s", "t"){ 
				#add seq
				$sf_stats{$subfam}{$s_or_t}{"seqs"}{$seqs->{$s_or_t}}=0; 
				#add sites
				# print $seqs->{$s_or_t} ."\n"; 
				(my $chr, my $start, my $end, my $strand) = $seqs->{$s_or_t} =~ /^([^:]+):(\d+)-(\d+)([+-])$/; 
				my @sites_ar = split(',', $sites->{$s_or_t});
				foreach my $pos(@sites_ar){
					my $coord = $chr . ":" . ($strand eq '+' ? ($start + $pos) : ($end - $pos + 1) ); #get correct position based on + or - 
					# print $coord ."\n"; 
					$sf_stats{$subfam}{$s_or_t}{"sites"}{$coord}=0; 
				}
			}
		}
	} while($subfam);
	close($cf_fh);
	#convert IDs of each subfamily to count of IDs of each subfamily
	foreach my $sf (keys %sf_stats){ #replace subfam hashes with the NUMBER of editing/edited 
		$sf_stats{$sf}{"pairs"} = keys %{$sf_stats{$sf}{"pairs"}}; 
		foreach my $s_or_t ("s", "t"){
			foreach my $s_or_s("seqs", "sites"){
				my $val = keys %{$sf_stats{$sf}{$s_or_t}{$s_or_s}}; 
				$sf_stats{$sf}{$s_or_t}{$s_or_s} = keys %{$sf_stats{$sf}{$s_or_t}{$s_or_s}};
			}
		}
	}
	return (\%sf_stats);
}

#Function: reads from a cluster-file file-handle and returns the (1) subfam, (2) source and (3) target in the next cluster
#			If no cluster was found (EOF reached) - returns NULL.
#Cluster tuple example: GA      anoCar2 LINE    L2      L2-1_ACar       chr3:75066880-75067118- chr1:49496-49711-       3       6       0       120,130,175     97,107,153      1       1,2,3   57      57      209     3|0|0|0|0|0|2|0|1|0|0|0
sub readTabCluster{
	my $fh = shift; 
	if (my $line = <$fh>){
		chomp $line; 
		my @f = split(/\t/, $line); 
		#get source, target and subfam from cluster
		my %seqs = ();
		my %sites = (); 
		(my $subfam, $seqs{"s"}, $seqs{"t"}, $sites{"s"}, $sites{"t"}) = ($f[4], $f[5], $f[6], $f[10], $f[11]); 
		return ($subfam, \%seqs, \%sites); 
	}
}

#prints two entries and the FP% (no endline)
#FP is the percent of the control of the real result
sub printWithFP{
	(my $fh, my $denom, my $nom) = @_; #file handle (open), denominator (real), nominator (control)
	print $fh "\t".($denom==0 ? '0' : $denom)."\t".($nom==0 ? '0' : $nom)."\t";
	my $FP = ($denom != 0 ) ? ($nom / $denom ) * 100 : 0;
	printf ($fh "%.3f", "$FP");
}

#Function: get Parameters from cluster_files_hash_ref
sub getParams{
	(my $cfh, my $paramLimits) = @_; 
	#get and reformat limits (reformat pval)
	(my $pval_h, my $pval_l, my $th_l, my $th_h) = split(',', $paramLimits); #8,16,3,16 for 8<=pvals<=16 and 3<=thresholds(min_clust_len)<=16  
	$pval_h = "1e-" . $pval_h if ($pval_h > 0 and $pval_h !~ /1e-/); 
	$pval_l = "1e-" . $pval_l if ($pval_l > 0 and $pval_l !~ /1e-/); 
	my %paramPairs = (); 
	foreach my $f (@{$cfh->{"GA"}}){ #Any parameter runned has a file in GA (and all mismatch-dirs)
		(my $org, my $class, my $fam, my $pval, my $th) = $f =~ /\S*clusters_([^_]+)_([^_]+)_(\S+)_([^_]+)_([^_]+).tab$/; 
		if ($paramLimits){
			next if ( ($pval_h and $pval > $pval_h) || ($pval_l and $pval < $pval_l) || ($pval_h and $pval > $pval_h) || ($pval_l and $pval < $pval_l) ); 
		}
		$paramPairs{$pval."\t".$th} = 0;
	}
	return \%paramPairs; 
}

