#Function: parses all cluster files of a family (of one org and class) and creates clust_stat files for families and subfamilies
#	Output: 1. Stat file for each family (regardless if any editing was found) - one file per family.
#			2. Stats for subfamilies that had at least some editing (Real or Control). 
#	Notes: 1. Removes all stat files currently existing (Family/Subfamily stats). 

use strict;
use getAll;  
(my $org , my $class, my $fam, my $subdir) = @ARGV; 
my $clustFiles = getAll::clustersFam($org, $class, $fam, 1, $subdir); 
die "no cluster_files for $org $class $fam\n" if $clustFiles == 0; 
(my $dir) = $clustFiles->[0] =~ /^(\S+)\/clusters_/; 
#remove old cluster_stats files
my $fam_for_cmd = $fam; $fam_for_cmd =~ s/\(/\\\(/g; $fam_for_cmd =~ s/\)/\\\)/g;
# system ("echo deleting clust_stats"); 
system("rm -f $dir/cluster_stats_".$org."_".$class."_".$fam_for_cmd.".txt");
# system ("echo deleting sf_clust_stats"); 
system("rm -f $dir/sf_cluster_stats_".$org."_".$class."_".$fam_for_cmd."_"."*".".txt");

my $subfams = getAll::names($org, $class, $fam); 

foreach my $cf (@$clustFiles){
	next if $cf =~ /control.txt$/;
	#extract data from file name
	(my $o, my $c, my $f, my $pval, my $th) = $cf =~ /clusters_([^_]+)_([^_]+)_(\S+)_(1e-\d+)_(\d+).txt$/;
	(my $editingPairs, my $edited) = getSubfamStats($cf); 
	# print (keys %$editingPairs) ."\t". (keys %$edited) ."\n";
	$cf =~ s/.txt$/_control.txt/; 
	(my $editingPairs_control, my $edited_control) = getSubfamStats($cf); 
	
	# Print subfamily cluster stats to subfamily stat file
	my $fam_editingPairs=0; my $fam_editingPairs_control=0; my $fam_edited=0; my $fam_edited_control=0;
		
	foreach my $sf(@$subfams){
		# Print subfamily cluster stats to subfamily stat file
		my $sf_stats = "$dir/sf_cluster_stats_".$org."_".$class."_".$fam."_".$sf.".txt";
		if (exists $editingPairs->{$sf} or exists $editingPairs_control->{$sf}){
			printStats($sf_stats, $org, $class, $fam, $pval, $th, $editingPairs->{$sf} , $editingPairs_control->{$sf} , $edited->{$sf} , $edited_control->{$sf}, $sf) ; 
		}
		# Sum data for family
		$fam_editingPairs += $editingPairs->{$sf}; 
		$fam_editingPairs_control += $editingPairs_control->{$sf}; 
		$fam_edited += $edited->{$sf}; 
		$fam_edited_control += $edited_control->{$sf}; 
	}
	#print stats for family
	my $fam_stats = "$dir/cluster_stats_".$org."_".$class."_".$fam.".txt"; 
	printStats($fam_stats, $org, $class, $fam, $pval, $th, $fam_editingPairs , $fam_editingPairs_control , $fam_edited , $fam_edited_control ); 
}

#Input: a clusters filename 
#Return: Pair of hashes: 	1. editingPairs{subfam} = num_unique_editing_pairs
#							2. edited{subfam} = num_unique_edited_sequences
#Note: The hashes returned contain only the subfams that had clusters 
sub getSubfamStats{
	(my $cf_name) = @_; 
	my %editingPairs = (); 
	my %edited = ();
	my $subfam, my $source, my $target; 
	open (my $cf_fh, "<". $cf_name) || die "didn't open $cf_name\n"; 
	do{
		($subfam, $source, $target) = readCluster($cf_fh);
		if ($subfam){
			$editingPairs{$subfam}{$source ."\t". $target}=0;
			$edited{$subfam}{$target}=0; 
		}
	}
	while($subfam);
	close($cf_fh); 
	#convert IDs of each subfamily to count of IDs of each subfamily
	foreach my $sf (keys %editingPairs){ #replace subfam hashes with the NUMBER of editing/edited 
		$editingPairs{$sf} = keys %{$editingPairs{$sf}};
		$edited{$sf} = keys %{$edited{$sf}}; 
	}
	return (\%editingPairs, \%edited);
}

#Function: reads from a cluster-file file-handle and returns the (1) subfam, (2) source and (3) target in the next cluster
#			If no cluster was found (EOF reached) - returns NULL.
sub readCluster{
	my $fh = shift; 
	my $line = '';
	my $cluster = ''; 
	my $source; 
	my $target; 
	#read cluster
	do{
		$line = <$fh>;
		$cluster .= $line; 
	}while ($line ne '' and $line !~ /^End cluster\n$/); 
	
	#get source, target and subfam from cluster
	if ($cluster ne ''){
		($source, $target) = $cluster =~ /\n[ACTG] name = (\S+)\n[ACTG] name = (\S+)/; 
		(my $subfam) = ($source =~ /\S+=(\S+)/); 
		return ($subfam, $source, $target); 
	}
	return; 
}

sub printStats{
	(my $stats_file, my $org, my $class, my $fam, my $pval, my $th, my $editing, my $editing_control, my $edited, my $edited_control, my $sf) = @_;
	#for output to be correct - convert empty strings to zeros
	$editing = 0 if $editing == 0; 
	$editing_control = 0 if $editing_control == 0; 
	$edited = 0 if $edited == 0; 
	$edited_control = 0 if $edited_control == 0; 
	#print output
	open(STATS , ">>".$stats_file) || die ("$stats_file didn't open\n");
	print STATS $org."\t".$class."\t".$fam."\t". ($sf ? $sf."\t" : '') .$pval."\t".$th."\t";
	print STATS $editing."\t".$editing_control."\t"; 
	my $FP_editing = ($editing != 0 ) ? ($editing_control / $editing ) * 100 : 0; 
	my $FP_edited = ($edited != 0 ) ? ($edited_control / $edited ) * 100 : 0;
	printf (STATS "%.3f", "$FP_editing"); 
	print STATS "\t".$edited."\t".$edited_control."\t";
	printf (STATS "%.3f", "$FP_edited");
	print STATS "\n"; 
	close(STATS);
}