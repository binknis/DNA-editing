#Function: get divergence from consensus scores (mismatches, ins, del and total) for each coord in a file
#Input: 1.
#			rmsk table fields: [0]class, [1]fam, [2]subfam, [3]coords(chrN:\d+-\d+[+-] format), [4]swScore, [5]milliDiv, [6]milliDel, [7]milliIns
#Output: stdout with fields: [0]assembly, [1]class, [2]fam, [3]subfam, [4]swScore, [5]pcDiv, [6]pcDel, [7]pcIns, [8]totalDiff, [9]num_sf_seqs
use strict; 

(my $coordsFile, my $rmskTableDir) = @ARGV; 
my %coords = (); 
my $colsPreScores = 4; 
my $scoreCols = 4; 
open (COORDS, $coordsFile) or die "open $coordsFile\n"; 
while(my $l = <COORDS>){
	chomp $l; 
	my @f = split(/\t/, $l); 
	$coords{$f[0]."\t".$f[1]} = 1; 
}
close(COORDS); 

opendir(DIR, $rmskTableDir) or die "opendir $rmskTableDir\n"; 
my @files =  readdir(DIR);  #grep {/^table_LTR.*/} readdir(DIR); 
closedir(DIR); 

foreach my $file (@files){
	my $rmskTableFile = $rmskTableDir ."/". $file; 
	(my $assembly) = $file =~ /table_LTR_([^_]+)\.rmsk\.txt/;
	if ($rmskTableFile =~ /\.gz$/){
		open (TABLE, "gunzip -c $rmskTableFile |") or die "open $rmskTableFile\n"; 
	}
	else{
		open (TABLE, $rmskTableFile) or die "open $rmskTableFile\n"; 
	}
	while (my $l = <TABLE>){
		chomp $l; 
		my @f = split(/\t/, $l); 
		if ($coords{$assembly."\t".$f[3]}){
			print $assembly ."\t". $l ."\n";
		}	
	}
	close(TABLE); 
}
