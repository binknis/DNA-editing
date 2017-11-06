#Computes cluster stats for all clusters in an Organism's Class
#cluster stats = File with graph statistics
#OUTPUT files are created in a subdirectory "Graph/*" in the original cluster (input) dir.
#DEFAULT: cluster input dir - "../Data/Organism/Class/results"; graph output dir: "../Data/Organism/Class/results/Graph"
# labelFlag := 0 - complete label; 1 - "id=name"; 2 - "name=coords"; 3 - "id=coords=name"; (coords = "chr:start-end[+-]")
#nonDefaultDir := an alternative subdirectory for the input and output. Built as such: "../Data/Organism/Class/results/nonDefaultDir"
#NOTE: no % is printed, even in FP % fields
#NOTE: families can't have "_" signs.
#IMPORTANT NOTE: sync any changes with specificGraphCreator.pl for consistency. 
use strict ;
use lib 'graph_scripts'; 
use GraphTools;
use File::Path qw(mkpath); 

(my $organism , my $class, my $dontMakeGraphFile, my $labelFlag, my $nonDefaultDir) = @ARGV;
my $pvalue;
my $th;
my $family;
my $graph_dir;
my $edgeCount;
my $editedCount; 
my $edgeCount_control;
my $editedCount_control; 
my $FP_percent_edges;
my $FP_percent_edited_seq;
my $name;
my $path_clusters;
my $path_clusters_control;

my $dir = "../Data/". $organism . "/" .$class ."/results";
if ($nonDefaultDir ne ""){
	$dir = $dir ."/". $nonDefaultDir;
}

#delete all old cluster-stats files
unlink(<$dir/cluster_stats*>);

opendir(CLUSTERS, $dir) || die "$dir didn't open\n";
my @clusterFileList = sort{lc($a) cmp lc($b)}(readdir(CLUSTERS));
shift(@clusterFileList) while ($clusterFileList[0] =~ /^\./); #erase "." and ".." links
closedir(CLUSTERS);

foreach my $clusterFile (@clusterFileList)
{
 if($clusterFile =~ /clusters_(\S+)/)
 {
	next if $clusterFile =~ /control.txt$/; 
    next if (-d $dir."/".$clusterFile); 
  #extract data from file name
  $clusterFile =~ /clusters_([^_]+)_([^_]+)_(\S+)_(1e-\d+)_(\d+).txt$/;
  $family = $3;
  $pvalue = $4;
  $th = $5;
  
  #create Graph file (real)
  # $graph_dir = "tracks_".$1."_".$2."_".$3;
  # $graph_dir = $dir."/".$graph_dir ;
  $graph_dir = $dir . "/Graphs";
  $path_clusters =  $dir."/".$clusterFile;
  mkpath($graph_dir) unless $dontMakeGraphFile;
  ($edgeCount, $editedCount) = GraphTools::clustToGraph($graph_dir, $path_clusters, $dontMakeGraphFile, $labelFlag);
  
  #create Graph file (real)
  $clusterFile =~ /(\S+)\.txt/; 
  my $prefix = $dir."/".$1;
  $path_clusters_control =  $prefix . "_control.txt";
  $graph_dir = $graph_dir . "_control";
  mkpath($graph_dir) unless $dontMakeGraphFile;
  ($edgeCount_control, $editedCount_control) = GraphTools::clustToGraph($graph_dir , $path_clusters_control, $dontMakeGraphFile, $labelFlag); 
  
  #calculate FP rates
  $FP_percent_edges = ($edgeCount != 0 ) ? ($edgeCount_control / $edgeCount ) * 100 : 0; 
  $FP_percent_edited_seq = ($editedCount != 0 ) ? ($editedCount_control / $editedCount ) * 100 : 0; 
  
  #create cluster-stats field
  open(my $stat_file , ">>$dir/cluster_stats_".$organism."_".$class."_".$family.".txt") || die ("the cluster_stat file is not open");
  print $stat_file $organism."\t".$class."\t".$family."\t".$pvalue."\t".$th."\t";
  print $stat_file $edgeCount."\t".$edgeCount_control."\t"; 
  printf ($stat_file "%.3f", "$FP_percent_edges"); 
  print $stat_file "\t".$editedCount."\t".$editedCount_control."\t"; 
  printf ($stat_file "%.3f", "$FP_percent_edited_seq"); 
  print $stat_file "\n"; 
  
  close($stat_file);  
 }
} 
