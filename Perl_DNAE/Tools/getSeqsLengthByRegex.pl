#function: extracts a list of fasta sequences from a Fasta file 
#input: 1. The file (full path). 2. a list of regexes
#output: prints the sequences to screen. 

use strict; 
# use Bio::SeqIO;

(my $fasta) = @ARGV; 
(my @regexes) = @ARGV; 

open (FASTA, $fasta) || die "couldn't open $fasta\n"; 
my $seq = ""; 
my $seqMatches=0; 
while (my $line = <FASTA>){
	if ($line =~ /^>/){
		chomp $line; 
		my $defline = $line; 
		$seqMatches=0; 
		foreach my $reg (@regexes){
			$seqMatches = 1 if $defline =~ /$reg/; 
		}
		$seq = "";
		while ($line = <FASTA>){
			if ($line =~ /^>/){
				seek(FASTA, -length($line), 1);
				last;
			}
			$seq .= $line; 
		}
		$defline =~ s/(>|\r)//g; 
		$seq =~ s/(\s+)//g;
		my $len = length($seq); 
		if ($seqMatches){
			print $defline . "\t" . $len ."\n"; 
		}
	}
}
close (FASTA);
