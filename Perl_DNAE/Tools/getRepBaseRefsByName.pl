#function: extracts a list of fasta sequences from a Fasta file 
#input: 1. The file (full path).  2. a list of Repbase IDs.
#output: prints the sequences to screen. 

use strict; 

my $fasta = shift; 
(my @names) = @ARGV; 
my %idHash = (); 
foreach my $name (@names){
	$nameHash{$name} = 1; 
}

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
	(my $seqName) = ($seq =~ /^>([^=]+)=/); #get seqs ID from head of defline
	chomp $seq; 
	print $seq . "\n" if exists $nameHash{$seqName}; 
}
close (FASTA); 
