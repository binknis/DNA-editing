#function: extracts a list of fasta sequences from a Fasta file 
#input: 1. The file (full path). 
#		2. a list of IDs for each sequence (can't contain whitespaces).
#output: prints the sequences to stdout. 

use strict; 

(my $fasta, my @regexes) = @ARGV; 
#create regex (by regex-flag)


#fetch sequences with matching IDs
open (FASTA, $fasta) || die "couldn't open $fasta\n"; 
my $seq; 
while (my $line = <FASTA>){
	if ($line =~ /^>/){
		$seq = $line; 
		while ($line = <FASTA>){
			if ($line =~ /^>/){
				seek(FASTA, -length($line), 1);
				last;
			}
			$seq .= $line; 
		}
	}
	(my $seqID) = ($seq =~ /^>(.+)?\n/); #get seqs ID from head of defline
	chomp $seq; 
	# print $seqID ."\n"; 
	#print the sequence only if it matches one of the regexes
	my $found = 0;
	foreach my $regex(@regexes){
		$found = 1 if  $seqID =~ /$regex/;
	}
	print $seq . "\n" unless $found;
}
close (FASTA); 
