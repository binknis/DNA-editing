#Function: for each organism, and classes LINE, LTR, SINE, DNA, RC - extracts the best parameters from cluster_stats
# 'best' = 1st tier: max num_edited within marginal diff from target value for FP ($percentFP). 
#		   2nd tier: closest to  FP edited (unique edited elements). 
#	by default has internal script searching for 10% and 0%. (To change, this script must be modified). 
#ARGS: control := finds 'best' for control. 

#OUTPUT: 

#NOTES: if you want to screen a length-filter other than 80 you need to change $percent_len asignment. 

use strict; 
use getAll; 
use File::Path qw(mkpath); 
use Math::BigInt;

my $margin = 1; 
getBestLinesAllOrgsAndClasses(0,10, $margin); 
getBestLinesAllOrgsAndClasses(0,0, $margin); 
getBestLinesAllOrgsAndClasses("G", 10, $margin); 
getBestLinesAllOrgsAndClasses("G", 0, $margin); 
getBestLinesAllOrgsAndClasses("GA", 10, $margin); 
getBestLinesAllOrgsAndClasses("GA", 0, $margin); 

#Finds the best params for each family in all orgs and classes
#searches first for the best within the marginal value; if not found- finds params closest to the target FP value. 
sub getBestLinesAllOrgsAndClasses{
	(my $GA, my $percentFP, my $margin, my $control) = @_; 
	my $percent_len = 80; #const of length cutoff used when file extraction was executed
	my $bestLine;
	my $resultsDir = "../DNA_editing_results/bestParams";
	my $subdir;
	if ($GA eq "G" ||  $GA eq "GA"){
		$subdir = "Clusters_".$percent_len."_percent_length_" . $GA; 
	}
	else{
		$GA = "raw";
		$subdir = 0; 
	}
	my $dir = $resultsDir ."/results_".$GA."/".$percentFP."FP_edited"; 
	mkpath($dir); 
	my $results = $dir."/results_all_".$GA."_".$percentFP."FP.txt";
	# $results =~ s/\.txt$/_control.txt/ if $control;  
	open (RESULTS, ">$results") || die "couldn't open $results\n";
	my $errorFile = $resultsDir. "/errors_gettingBestLines.txt"; 
	open (ERRORS, ">>$errorFile") || die "couldn't open $errorFile\n";
	
	my $orgs = getAll::organisms(); 
	foreach my $org (@$orgs){
		my $classes = getAll::classes($org); 
		foreach my $class (@$classes){
			next unless $class =~ /^(LINE|LTR|SINE|DNA|RC|Other)$/; #***
			my $fams = getAll::families($org, $class); 
			foreach my $family (@$fams){
				my $clust_stats = getAll::clustStatsFam($org, $class, $family, 1, $subdir);
				unless ($clust_stats){ #no clustStats for this family (in this organism's class)
					print ERRORS "$org $class $family doesn't have cluster_stats in $subdir\n"; 
					next; 
				}
				foreach my $clust_stats_file (@$clust_stats){
					($bestLine) = getBestLine($clust_stats_file, $percentFP, $margin); 
					if ($bestLine eq ""){ #no lines have FP within margin - find closest which is outside of margin.
						($bestLine) = getBestLineNoMargin($clust_stats_file, $percentFP, $margin);
					}
					chomp $bestLine; 
					print RESULTS $org."\t".$class."\t".$bestLine ."\t". $percentFP ."\t". $GA . "\n";
				}
			}
		}
	}
	close(ERRORS);
	close(RESULTS); 
}

#finds line whose params produced highest num_edited and a FP within marginal difference from target FP value
sub getBestLine{
	(my $clust_stats_file, my $percentFP, my $margin) = @_; 
	my $max_edited = -1; 
	my $bestLine; 
	my $smallestDelta = Math::BigInt->binf(); 
	#get lines of cluster_stats file
	my $clust_stats = getAll::lines($clust_stats_file); 
	return 0 unless $clust_stats; #error
	#find line closest to the wanted percentage
	foreach my $line (@$clust_stats){
		(my $family, my $pval, my $th, my $edges, my $edges_control, my $fp_edges, my $edited, my $edited_control, my $fp_edited) = split(/\s+/,$line); 
		#if 'edited' is closer to the target percent, or equal and has higher rates of editing - save it as best
		if(abs($fp_edited-$percentFP) < $margin && $edited > $max_edited){
			$bestLine = $line; 
			$max_edited = $edited; 
			$smallestDelta = abs($fp_edited-$percentFP); 
		}
	}
	return $bestLine; 
}

#finds line with closest FP to target FP value (and by max num_edited). Doesn't filter by margin at all. 
sub getBestLineNoMargin{
	(my $clust_stats_file, my $percentFP) = @_; 
	my $max_edited = -1; 
	my $bestLine; 
	my $smallestDelta = Math::BigInt->binf(); 
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