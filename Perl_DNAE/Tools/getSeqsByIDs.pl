#function: extracts a list of fasta sequences from a Fasta file 
#input: 1. The file (full path). 
#		2. Regex flag (num= numerical id, def= whole defline(no white-spaces)
#		3. a list of IDs for each sequence (can't contain whitespaces).
#output: prints the sequences to stdout. 

use strict; 

(my $fasta, my $regexFlag , my @ids) = @ARGV; 
my %idHash = (); 
foreach my $id (@ids){
	$idHash{$id} = 1;
}

#create regex (by regex-flag)
my $regex;
if ($regexFlag eq "num"){
	$regex = '^>(\d+)';
}
elsif($regexFlag eq "def"){
	$regex = '^>(\S+)';
}
else{
	die "bad regex flag for getSeqByIDs.pl! must be num or def\n"; 
}

#fetch sequences with matching IDs
open (FASTA, "<$fasta") || die "couldn't open $fasta\n"; 
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
	(my $seqID) = ($seq =~ /$regex/); #get seqs ID from head of defline
	# (my $seqID) = ($seq =~ />(\d+)/);
	chomp $seq; 
	# print $seqID ."\n"; 
	print $seq . "\n" if exists $idHash{$seqID};
}
close (FASTA); 
