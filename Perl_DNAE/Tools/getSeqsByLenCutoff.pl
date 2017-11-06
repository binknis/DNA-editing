#function: extracts a list of fasta sequences from a Fasta file 
#input: 1. The file (full path). 2. a list of numerical IDs for each sequence.
#output: prints the sequences to screen. 

use strict; 
use Bio::SeqIO;

(my $fasta, my $lenCutoff) = @ARGV; 

open (FASTA, $fasta) || die "couldn't open $fasta\n"; 
my $seq = ""; 
my $seqID; 
while (my $line = <FASTA>){
	if ($line =~ /^>/){
		($seqID) = ($line =~ /^>(.*)/); #get seqs ID from head of defline 
		$seq = ""; 
		while ($line = <FASTA>){
			if ($line =~ /^>/){
				seek(FASTA, -length($line), 1);
				last;
			}
			$seq .= $line; 
		}
		$seq =~ s/(\s+)//g;
		if(length($seq) >= $lenCutoff){
			print ">" . $seqID ."\n". $seq ."\n"; 
		}
	}
}
close (FASTA);
