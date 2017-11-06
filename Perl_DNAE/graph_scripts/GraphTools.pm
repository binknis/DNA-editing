package GraphTools; 
use strict;

### subroutine getNodeLabel (used by sub "clustToGraph") ###
#return node label by specific parameters
# flag := 0 - complete label; 1 - "id=name"; 2 - "name=coords"; 3 - "id=coords=name"; (coords = "chr:start-end[+-]")
#sequence description format example: 2410=panTro3=chr20:47198565-47200078+=Other=Other=SVA_D
sub getNodeLabel{
	my $label; 
	(my $line, my $flag) = @_;
	$line =~ /[AG] name = (\S+)/;
	$line = $1; 
	my @fields = split(/=/, $line); 
	if ($flag == 0){ #label is full label
	   $label = $line;
	}
	elsif ($flag == 1){ #label is id=name (together are unique)
	   $label = $fields[0] . "=" . $fields[5]; 
	}
	elsif ($flag == 2){ #label is coordinates=name
		$label = $fields[2] . "=" . $fields[5];
	}
	elsif ($flag == 3){ #label is id=coords=name
		$label = $fields[0] . "=" . $fields[2] . "=" . $fields[5]; 
	}
	else {
		print "you have a problem with your label-flag!\n"; 
		exit; 
	}
	return $label; 
}

### subroutine clustToGraph ###
#Creates a graph file from a cluster file and returns 
#INPUT: cluster file.
#OUTPUT: list of the sources and targets to "graph.txt" file,
#in addition, printing of amount of uniqe sources and targets to screen.
#dontMakeFile this is flag if to create graph file or not.

#Arguments: 
#graph_dir = output dir for graph
#clusterFileFullPath = FULL PATH of cluster file, including file name.
# dontMakeFile = a flag. 0: make file. 1: don't. 
#labelFlag = determines which fields to use for annotation in graph file (see sub "getLabel")

#Notes:
# cluster file format: "G name = ... " and in next line "A name = ... "
sub clustToGraph 
{
	(my $graph_dir, my $clusterFileFullPath, my $dontMakeFile, my $labelFlag) = @_; 
	open (CLUSTER_FILE ,$clusterFileFullPath) || die ("couldn't open cluster file");  

	my $graphFile; my $graphFileName;
	unless ($dontMakeFile){
		$clusterFileFullPath =~ /clusters_([^_]+)_([^_]+)_(\S+)_(1e-\d+)_(\d+)(_control)?\.txt$/;
		$graphFileName = "$graph_dir/graph_$1_$2_$3_$4_$5$6.txt";
		open ($graphFile, ">$graphFileName")  || die ("couldn't open graph file");
	}

	my %edges = (); #will hold all unique edges (G->A). (Needed because every HSP generates a cluster, even if editing was found in two HSPs of a single result). 
	my %editedNodes = (); #will hold all unique targets (A)
	my $source;  #G
	my $target;  #A
	my $g;
	my $a;
	my @arrayG = ();
	my @arrayA = ();

	#read cluster file and insert edges and edited nodes into their hashes. 
	while(my $line  = <CLUSTER_FILE>)
	{
	  chomp;
	  if($line =~/^Found/)
		{
			#save G-name
		  $g  = <CLUSTER_FILE>; 
		  $source = &getNodeLabel($g, $labelFlag); 
		  push(@arrayG,$source);		
		  #get A-name
		  $a  = <CLUSTER_FILE>; 
		  $target = &getNodeLabel($a, $labelFlag); 
		  $editedNodes{$target} = 1; 
		  push(@arrayA,$target);
		  #print to file
		  print $graphFile "$source\t$target\n" unless ($dontMakeFile || exists $edges{$source."\t".$target});		  
		  $edges{$source."\t".$target} = 1;
		}
	}
	close(CLUSTER_FILE);
	close($graphFile) unless $dontMakeFile; 
	
	my $edgeCount  = scalar keys (%edges) ;  
	my $editedCount = scalar keys (%editedNodes); 
	return ($edgeCount, $editedCount);
}



1;