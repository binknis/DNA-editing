#Function: retains a specific list of coord-pairs based on coord file
#SorT = where to apply filter: 0 - source (1st column), 1 - target (2nd column)
#Input seq file could be with or without strand specifier
use strict; 

(my $pairFile, my $seqFile, my $SorT) = @ARGV; 

my %seqs = (); 
my $noStrand = 0; 
open (SEQS ,$seqFile) || die ("couldn't open $seqFile"); 
while (my $l = <SEQS>){
	chomp $l; 
	my $coords; 
	if ($l =~ /=/){
		($coords) = $l =~ /=([^=]+:\d+-\d+[+-])/;
		$coords =~ s/\s+//g; 
	}
	else{
		$coords = $l;
	}
	$seqs{$coords}=0;
	$noStrand = 1 if $coords !~ /\S+:\d+-\d+[+-]/;
}
close(SEQS); 

open (PAIRS ,$pairFile) || die ("couldn't open $pairFile"); 
while (my $l = <PAIRS>){
	chomp $l; 
	$l =~ s/=\S+//g; 
	my @pair = split (/\t/, $l); 
	if ($noStrand){
		$pair[0] =~ s/[+-]$//; 
		$pair[1] =~ s/[+-]$//; 
	}
	next if ($SorT==0 and not exists $seqs{$pair[0]}); #filter on source
	next if ($SorT==1 and not exists $seqs{$pair[1]}); #filter on target
	print $l ."\n"; 
}
close(PAIRS);
