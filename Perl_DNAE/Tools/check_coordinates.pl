#check if the coordinates in a interval file are equal to the coords in a sequence file
#prints the amount of coords in each file that weren't found in it's counterpart

use strict; 
( my $seqDataFile, my $seqFile) = @ARGV;

my %coords = (); 
open(SEQDATA,"<$seqDataFile") || die ("couldn't open description file\n"); 
my $chrom; 
my $start; 
my $end; 
my $description; 
while (my $line = <SEQDATA>){
	chomp; 
	if ($line =~ /^\s*#/){next;} #skip remark lines
	my @data = split(/\s+/,$line); 
	$chrom = $data[0]; 
	$start = $data[1]; 
	$end = $data[2]; 
	$description = $chrom."_".$start."_".$end; 
	unless (exists $coords{$description}){
		$coords{$description}=0; 
	}
	$coords{$description}++; 
}
close(SEQDATA); 

open(SEQS,"<$seqFile") || die ("couldn't open sequence file\n"); 
while (my $line = <SEQS>){
	if ($line !~ "^>"){next;}
	chomp;
	$line =~ (/^>\S+_(chr\S+)_(\d+)_(\d+)_[+-]\s*$/);
	$chrom = $1; 
	$start = $2; 
	$end = $3; 
	$description = $chrom."_".$start."_".$end; 
	unless (exists $coords{$description}){
		$coords{$description}=0; 
	}
	$coords{$description}--; 
}
close(SEQS); 

#print results to files
open (my $extraData, ">extraData.txt") || die "couldn't open extraData.txt"; 
open (my $extraSeq, ">extraSeq.txt") || die "couldn't open extraSeq.txt"; 
#sum amount of coords > 0 and coords < 0 seperately
 while ( my ($coord, $num) = each(%coords) ) {
	if ($num > 0){
		print $extraData "$coord: $num\n";
	}
	elsif ($num < 0){
		print $extraSeq "$coord: $num\n";
	}
 }
close ($extraData);
close ($extraSeq);




