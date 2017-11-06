use strict; 
use File::Path qw(mkpath); 
use getAll; 

(my $class, my $family, my $percentFP, my $subdir) = @ARGV; 
my $best_pval; 
my $best_th; 
my $graphDir = ""; 
mkpath($graphDir);

my $dataDir = "../Data"; 
my $orgs = getAll::organisms();

foreach my $org (@$orgs){
	my $resultDir = "$dataDir/$org/$class/results/$subdir";
	my $maxEdited = -1;
	$best_pval = -1;
	$best_th = -1;
	#read organism's clustStats file
	my $clustStats = $resultDir . "/cluster_stats_".$family.".txt";
	next unless (-d $resultDir && -f $clustStats); #skip organisms that don't have the subdirectory (may not have LINE at all)
	print "$org\n"; #***
	open(STATS, $clustStats) || print "$clustStats didn't open\n"; 
	while (my $line = <STATS>){
		(my $family, my $pval, my $th, my $edges, my $edges_control, my $fp_edges, my $edited, my $edited_control, my $fp_edited) = split(/\s+/,$line); 
		if($fp_edited <= $percentFP && $edited > $maxEdited){
			$maxEdited = $edited; 
			$best_pval = $pval; 
			$best_th = $th; 
		}
	}
	close(STATS); 
	
	### create graph file ###
	#cluster file format: clusters_Cat_LINE_CR1_1e-10_10_control.txt
	#(my $clusterFullPath, my $outDir, my $dontMakeFile, my $labelFlag, my $printOutput) = @ARGV; 
	my $clustFile = $resultDir."/clusters_".$org."_".$class."_".$family."_".$best_pval."_".$best_th.".txt";
	system("perl specificGraphCreator.pl $clustFile $graphDir 0 3 &");
	# my $clustFile_control = $resultDir."/clusters_".$org."_".$class."_".$family."_".$pval."_".$th."_control.txt";
	# system("specificGraphCreator.pl $clustFile_control $graphDir 0 3 &");	
}