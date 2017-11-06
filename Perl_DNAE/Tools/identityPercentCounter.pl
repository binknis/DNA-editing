#Input: Blast output file. 
#Function: prints histograms of - 1. The PERCENTAGE of identity in every alignment. 
# 								  2. The AMOUNT of mismatches in every alignment. 
#extracts all data from "Identities" line. example: "Identities = 148/148 (100%)"

use strict; 

my $blastFile = shift; 
open(BLAST, $blastFile) || die "didn't open $blastFile\n"; 
my %identity = (); 
my %unidentical = (); 
my $diff; 
my $totalAlignments = 0; 

while (my $line = <BLAST>){

	if ($line =~ /^\s*Identities = (\d+)\/(\d+) \((\d+)%\)/){
		$totalAlignments++; 
		#save amount of unidentical nucs
		$diff = $2 - $1; 
		unless (exists $unidentical{$diff}){
			$unidentical{$diff} = 0;
		}
		$unidentical{$diff}++; 
		#save percentage of identity
		unless (exists $identity{$3}){
			$identity{$3} = 0; 
		}
		$identity{$3}++;
	}
}
#print total amount of alignments
print "Total amount of alignments: $total\n"; 

#print unidentical nuc count
foreach my $key(sort{$b <=> $a}(keys %unidentical)){
	print $key . "\t" . $unidentical{$key} . "\n"; 
}

#print identity percent count
my $over90Count = 0; 
foreach my $key(sort{$b <=> $a}(keys %identity)){
	print $key . "\t" . $identity{$key} . "\n"; 
	if ($key > 90){
		$over90Count += $identity{$key}; 
	}
}
my $over90Percent = $over90Count / $total * 100; 
print "Percent of alignments over 90% identity: $over90Percent\n";

