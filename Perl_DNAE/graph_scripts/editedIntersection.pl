#INPUT: two cluster files
#the script check if exist same name of cluster in cluster_file and cluster_file_control.
#OUTPUT: Two files containing the As and Gs found in both files.
#NOTE: I didn't check that it works. 
use strict;
(my $clusterFile1, my $clusterFile2, my $outFilePrefix) = @ARGV;
my %commonA = ();
my %commonG = ();
my $a;
my $g;
my $commonAFile; 
my $commonGFile; 

open (CLUST1 ,$clusterFile1) || die ("couldn't open cluster file"); 
if ($outFilePrefix){
	$clusterFile1 =~ /_([^_]+)_(1e-\d+)_(\d+).txt$/;
	open ($commonAFile, ">intersectA_".$1."_".$2."_".$3.".txt");
	open ($commonGFile, ">intersectG_".$1."_".$2."_".$3.".txt");
}
else{
	open ($commonAFile, ">$outFilePrefix"."_A.txt");
	open ($commonGFile, ">$outFilePrefix"."_G.txt");
}

while(my $line  = <CLUST1>)
{
  chomp;
  if($line =~/^Found/)
	{
	  chomp($a  = <CLUST1>);			 
	  $commonA{$1.$2} = 0;
	  chomp($g  = <CLUST1>);		
	  $commonG{$1.$2} = 0;		 		
	}  
}	
close(CLUST1);

open (CLUST2 ,$clusterFile2) || die ("couldn't open cluster file"); 
while(my $line  = <CLUST2>)
{
  chomp;
  if($line =~/^Found/)
	{
	  chomp($a  = <CLUST2>);			
	  $a =~ /A name = (\d+=).+=(\S+)/;
	  if (exists  $commonA{$1.$2}){ 
		 $commonA{$1.$2}++;
	  }
	  
	  chomp($g  = <CLUST2>);
	  $g =~ /G name = (\d+=).+=(\S+)/;
	  if (exists  $commonG{$1.$2}){ 
		$commonG{$1.$2}++;
	  }		  		 		
	} 
}	
close(CLUST2);  


#print intersecting As to file
foreach $a (keys(%commonA)){
	print $commonAFile "$a\n";
}
foreach $g (keys(%commonG)){
	print $commonGFile "$g\n";
}
