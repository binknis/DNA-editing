#Function: parses all cluster files of a family (of one org and class) and creates clust_stat files for families and subfamilies
#	Output: 1. Stat file for each family (regardless if any editing was found) - one file per family.
#			2. Stats for subfamilies that had at least some editing (Real or Control). 
#		Output fields: Org	class	fam	subfam	pval th	edit_pairs	edit_pairs_c	FP_pairs	edited	edited_c	FP_edited	
#							----> edited_sites	edited_sites_c	FP_edited_sites	source	source_c	FP_source source_sites	source_sites_c FP_source_sites
#	Notes: 1. Removes all stat files currently existing (Family/Subfamily stats). 
#started changing. Find a way to save sites{G} and sites{A}
#			2. Names must have 0-base starts.  

use strict;
use lib $ENV{"HOME"} ."/Perl_DNAE"; 
use getAll;  
(my $org , my $class, my $fam, my $subdir) = @ARGV; 
my $clustFiles = getAll::clustersFam($org, $class, $fam, 1, $subdir); 
die "no cluster_files for $org $class $fam\n" if $clustFiles == 0; 
(my $dir) = $clustFiles->[0] =~ /^(\S+)\/clusters_/; 
#remove old cluster_stats files
my $fam_for_cmd = $fam; $fam_for_cmd =~ s/\(/\\\(/g; $fam_for_cmd =~ s/\)/\\\)/g;
# system ("echo deleting clust_stats2"); 
system("rm -f $dir/cluster_stats2_".$org."_".$class."_".$fam_for_cmd.".txt");
# system ("echo deleting sf_clust_stats2"); 
system("rm -f $dir/sf_cluster_stats2_".$org."_".$class."_".$fam_for_cmd."_"."*".".txt");

my %sf_stats = (); #{real}, {ctrl}
my %fam_stats = (); 
my $subfams = getAll::names($org, $class, $fam); 

foreach my $cf (@$clustFiles){
	next if $cf =~ /control.txt$/;
	#extract data from file name
	(my $o, my $c, my $f, my $pval, my $th) = $cf =~ /clusters_([^_]+)_([^_]+)_(\S+)_(1e-\d+)_(\d+).txt$/;
	# (my $editingPairs, my $edited, my $editedSites) = getSubfamStats($cf); 
	$sf_stats{"real"} = getSubfamStats($cf); 
	
	# print (keys %$editingPairs) ."\t". (keys %$edited) ."\n";
	$cf =~ s/.txt$/_control.txt/; 
	# (my $editingPairs_c, my $edited_c, my $editedSites_c) = getSubfamStats($cf); 
	$sf_stats{"ctrl"} = getSubfamStats($cf); 
	
	#init %fam_stats HoHoHoH
	foreach my $r_or_c ("real", "ctrl") {
		$fam_stats{$r_or_c}{"pairs"} = 0; 
		foreach my $s_or_t ("s", "t"){
			foreach my $s_or_s("seqs", "sites"){
				$fam_stats{$r_or_c}{$s_or_t}{$s_or_s} = 0;
			}
		}
	}
	
	# Print subfamily cluster stats to subfamily stat file	
	foreach my $sf(@$subfams){
		# Print subfamily cluster stats to subfamily stat file
		my $sf_stats_file = "$dir/sf_cluster_stats2_".$org."_".$class."_".$fam."_".$sf.".txt";
		if (exists $sf_stats{"real"}{$sf} or exists $sf_stats{"ctrl"}{$sf}){
			printStats($sf_stats_file, $org, $class, $fam, $pval, $th, \%sf_stats, $sf) ;
		}
		# Sum data for family
		foreach my $r_or_c ("real", "ctrl") {
			$fam_stats{$r_or_c}{"pairs"} += $sf_stats{$r_or_c}{$sf}{"pairs"}; 
			foreach my $s_or_t ("s", "t"){
				foreach my $s_or_s("seqs", "sites"){
					$fam_stats{$r_or_c}{$s_or_t}{$s_or_s} += $sf_stats{$r_or_c}{$sf}{$s_or_t}{$s_or_s};
				}
			}
		}
	}
	#print stats for family
	my $fam_stats_file = "$dir/cluster_stats2_".$org."_".$class."_".$fam.".txt"; 
	printStats($fam_stats_file, $org, $class, $fam, $pval, $th, \%fam_stats); 
}

#Input: a clusters filename 
#Return: Hash sf_stats holds (for each subfam): a. pairs; b. source seqs; c. source sites; d. target seqs ; e. target sites. 
#Note: The hashes returned contain only the subfams that had clusters 
sub getSubfamStats{
	(my $cf_name) = @_; 
	my %sf_stats = (); #"pairs", {"s"}({seqs}, {sites}), {"t"}({seqs}, {sites}),
	my $subfam, my $seqs; my $sites; 
	open (my $cf_fh, "<". $cf_name) || die "didn't open $cf_name\n"; 
	do{
		($subfam, $seqs, $sites) = readCluster($cf_fh);
		if ($subfam){ #cluster fetched
			#add pair
			$sf_stats{$subfam}{"pairs"}{$seqs->{"s"} ."\t". $seqs->{"t"}}=0; 
			#Add seq and sites for source and target
			foreach my $s_or_t("s", "t"){ 
				#add seq
				$sf_stats{$subfam}{$s_or_t}{"seqs"}{$seqs->{$s_or_t}}=0; 
				#add sites
				# print $seqs->{$s_or_t} ."\n"; 
				(my $chr, my $start, my $end, my $strand) = $seqs->{$s_or_t} =~ /\d+=[^=]+=([^:]+):(\d+)-(\d+)([+-])=[^=]+=[^=]+=\S+/; 
				my @sites_ar = split(/\s+/, $sites->{$s_or_t});
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
				#print $sf."\t".$s_or_s."\t".$s_or_t."\t".$s_or_s."\t".(scalar keys %{$sf_stats{$sf}{$s_or_t}{$s_or_s}}) ."\n"; 
				# print $sf."\t".$s_or_t."\t".$s_or_s."\t".(scalar $val) ."\n"; 
			}
		}
	}
	return (\%sf_stats);
}

#Function: reads from a cluster-file file-handle and returns the (1) subfam, (2) source and (3) target in the next cluster
#			If no cluster was found (EOF reached) - returns NULL.
sub readCluster{
	my $fh = shift; 
	my $line = '';
	my $cluster = ''; 
	my $sites_s; 
	my $sites_t;
	#read cluster
	do{
		$line = <$fh>;
		$cluster .= $line; 
	}while ($line ne '' and $line !~ /^End cluster\n$/); 
	
	#get source, target and subfam from cluster
	my %seqs = ();
	my %sites = (); 
	if ($cluster ne ''){
		($seqs{"s"}, $seqs{"t"}) = $cluster =~ /\n[ACTG] name = (\S+)\n[ACTG] name = (\S+)/; 
		($sites{"s"}, my $ignore, $sites{"t"}) = $cluster =~ /Locations [ACTG] seq:\s+((\d+|\s+)+)\nLocations [ACTG] seq:\s+((\d+|\s+)+)\n/;
		 #print "source: ". $sites{"s"} ."\t". "target: ". $sites{"t"} . "\n"; 
		(my $subfam) = ($seqs{"t"} =~ /\S+=(\S+)/);
		 #print $seqs{"s"}."\t". $seqs{"t"}."\t". $sites{"s"}."\t".  $sites{"t"} ."\t". $subfam."\n";  #***
		return ($subfam, \%seqs, \%sites); 
	}
	return; 
}

#print output (for family or subfamily; depends on $sf param)
sub printStats{
	(my $stats_file, my $org, my $class, my $fam, my $pval, my $th, my $stats, my $sf) = @_;
	my $r; #real 
	my $c; #ctrl
	if ($sf){ #sf_stats
		$r = $stats->{"real"}{$sf}; 
		$c = $stats->{"ctrl"}{$sf}; 
	}
	else{ #fam_stats
		$r = $stats->{"real"}; 
		$c = $stats->{"ctrl"}; 
	}
	#print stats
	open(my $STATS , ">>".$stats_file) || die ("$stats_file didn't open\n");
	print $STATS $org."\t".$class."\t".$fam."\t". ($sf ? $sf."\t" : '') .$pval."\t".$th;
	printWithFP($STATS, $r->{"pairs"},$c->{"pairs"});
	
	foreach my $s_or_t ("t", "s"){ #"t" vs "s" are reversed in comparison to other loops
		foreach my $s_or_s("seqs", "sites"){ 
			printWithFP($STATS, $r->{$s_or_t}{$s_or_s},$c->{$s_or_t}{$s_or_s});
		}
	}
	print $STATS "\n";
	close($STATS);
}

#prints two entries and the FP% (no endline)
#FP is the percent of the control of the real result
sub printWithFP{
	(my $fh, my $denom, my $nom) = @_; #file handle (open), denominator (real), nominator (control)
	 
	print $fh "\t".($denom==0 ? '0' : $denom)."\t".($nom==0 ? '0' : $nom)."\t";
	my $FP = ($denom != 0 ) ? ($nom / $denom ) * 100 : 0;
	printf ($fh "%.3f", "$FP");
}