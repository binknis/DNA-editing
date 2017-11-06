use strict;
use lib 'graph_scripts'; 
use GraphTools; 

### subroutine clustToGraph ###
#Creates a graph file from a cluster file and returns 
#INPUT: cluster file.
#OUTPUT: creates a graph file and/or prints the amount of unique edges and edited sequences.
#dontMakeFile this is flag if to create graph file or not.
#NOTE: parsing is problematic for Organism, Class or Family with underscores

#example: perl graph_scripts/specificGraphCreator.pl ../Data/AnolisL2/LINE/results/clusters_AnolisL2_LINE_L2_1e-8_6.txt "../DNA_editing_results" 0 0 0 & 


(my $clustersFullPath, my $graph_dir, my $dontMakeFile, my $labelFlag, my $printOutput) = @ARGV; 

(my $edgeCount, my $editedCount) = GraphTools::clustToGraph($graph_dir, $clustersFullPath, $dontMakeFile, $labelFlag);

print ($edgeCount, $editedCount) if $printOutput;