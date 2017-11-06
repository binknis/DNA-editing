package GraphTools; 
use strict;

### subroutine clustToGraph ###
#Creates a graph file from a cluster file and returns 
#INPUT: cluster file.
#OUTPUT: list of the sources and targets to $graphFile file,
#in addition, printing of amount of uniqe sources and targets to screen.
#dontMakeFile this is flag if to create graph file or not.
#*** parsing is problematic for Organism, Class or Family with underscores


(my $clusterFullPath, my $outDir, my $dontMakeFile, my $printOutput) = @ARGV; 
my %edges = (); #will hold all unique edges (G->A). (Needed because every HSP generates a cluster, even if editing was found in two HSPs in of result). 
my %editedNodes = (); #will hold all unique targets (A)
my $target;  #A
my $source;  #G
my $a;
my $g;
my @arrayA = ();
my @arrayG = ();

$clusterFullPath =~ /\/([^\/]+)$/; 
my $clusterFileName = $1; 
$clusterFileName =~ /_([^_]+)_([^_]+)_([^_]+)_(1e-\d+)_(\d+)(_control)?\.txt$/; 
#$clusterFileName =~ /_([^_]+)_(1e-\d+)_(\d+)(_control)?\.txt$/; 
my $graphFile = "graph_$1_$2_$3_$4_$5$6.txt";

open (GRAPH, ">$outDir/$graphFile")  || die ("couldn't open $outDir/$graphFile");
open (CLUSTERS ,$clusterFullPath) || die ("couldn't open $clusterFullPath");
#read cluster file and insert edges and edited nodes into their hashes. 
while(my $line  = <CLUSTERS>)
{
  chomp;
  if($line =~/^Found/)
	{
	  $a  = <CLUSTERS>; 
	  $a =~ /A name = (\d+=).+=(\S+)/;   #get id and name (together are unique)
	  $target = $1.$2;		  		  
	  $editedNodes{$target} = 1; 
	  push(@arrayA,$target);
	  $g  = <CLUSTERS>; 
	  $g =~ /G name = (\d+=).+=(\S+)/;  #get id and name (together are unique)
	  $source = $1.$2;
	  push(@arrayG,$source);		
	  print GRAPH "$source\t$target\n" unless ($dontMakeFile || exists $edges{$source."\t".$target}); 
	  $edges{$source."\t".$target} = 1;
	}  
}

close(CLUSTERS);
close(GRAPH) unless $dontMakeFile; 

my $edgeCount  = scalar keys (%edges) ;  
my $editedCount = scalar keys (%editedNodes); 
print ($edgeCount, $editedCount) if $printOutput;
