#Function: for each organism, and classes LINE, LTR, SINE - extracts the best parameters from cluster_stats
#			file and creates the corellating graphs. 
# 'best' = closest to 10% and 0% FP edited (unique edited elements). 

#ARGS: control := finds 'best' for control. 

#OUTPUT: 
#	1. creates two graphs for each family

#NOTES: if you want to screen a length-filter other than 80 you need to change $percent_len asignment. 

use strict; 
use getAll; 
use File::Path qw(mkpath); 

getBestLinesAllOrgsAndClasses(0,10); 
getBestLinesAllOrgsAndClasses(0,0); 
getBestLinesAllOrgsAndClasses("G", 10); 
getBestLinesAllOrgsAndClasses("G", 0); 
getBestLinesAllOrgsAndClasses("GA", 10); 
getBestLinesAllOrgsAndClasses("GA", 0); 


sub getBestLinesAllOrgsAndClasses{
	(my $GA, my $percentFP, my $control) = @_; 
	my $percent_len = 80; #const of length cutoff used when file extraction was executed
	my $bestLine;
	my $resultsDir = "../DNA_editing_results/bestParams";
	my $subdir;
	if ($GA eq "G" ||  $GA eq "GA"){
		$subdir = "Clusters_".$percent_len."_percent_length_" . $GA; 
	}
	else{
		$GA = "rawClusters";
		$subdir = 0; 
	}
	my $dir = $resultsDir ."/results_".$GA."/".$percentFP."FP_edited"; 
	mkpath($dir); 
	my $results = $resultsDir."/results_all_".$GA."_".$percentFP."FP.txt";
	# $results =~ s/\.txt$/_control.txt/ if $control;  
	open (RESULTS, ">$results") || die "couldn't open $results\n";
	my $errorFile = $resultsDir. "/errors_gettingBestLines.txt"; 
	open (ERRORS, ">>$errorFile") || die "couldn't open $errorFile\n";
	
		
	my $orgs = getAll::organisms(); 
	foreach my $org (@$orgs){
		my $classes = getAll::classes($org); 
		foreach my $class (@$classes){
			next unless $class =~ /^(LINE|LTR|SINE)$/; #***
			my $fams = getAll::families($org, $class); 
			foreach my $family (@$fams){
				my $clust_stats = getAll::clustStatsFam($org, $class, $family, 1, $subdir);
				unless ($clust_stats){ #no clustStats for this family (in this organism's class)
					print ERRORS "$org $class $family doesn't have cluster_stats in $subdir\n"; 
					next; 
				}
				foreach my $clust_stats_file (@$clust_stats){
					($bestLine) = getBestLine($clust_stats_file, $percentFP); 
					chomp $bestLine; 
					print RESULTS $org."\t".$class."\t".$bestLine ."\t". $percentFP ."\t". $GA . "\n";
				}
			}
		}
	}
	close(ERRORS);
	close(RESULTS); 
}

sub getBestLine{
	(my $clust_stats_file, my $percentFP) = @_; 
	my $max_edited = -1; 
	my $bestLine; 
	my $smallestDelta = 100000000000000; 
	#get lines of cluster_stats file
	my $clust_stats = getAll::lines($clust_stats_file); 
	return 0 unless $clust_stats; #error
	#find line closest to the wanted percentage
	foreach my $line (@$clust_stats){
		(my $family, my $pval, my $th, my $edges, my $edges_control, my $fp_edges, my $edited, my $edited_control, my $fp_edited) = split(/\s+/,$line); 
		#if 'edited' is closer to the target percent, or equal and has higher rates of editing - save it as best
		if(abs($fp_edited-$percentFP) < $smallestDelta || (abs($fp_edited-$percentFP) == $smallestDelta && $edited > $max_edited)){
			$bestLine = $line; 
			$max_edited = $edited; 
			$smallestDelta = abs($fp_edited-$percentFP); 
		}
	}
	return $bestLine; 
}
