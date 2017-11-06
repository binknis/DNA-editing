#Function: Converts sequences in seqFile to FASTA format (number of characters per line can be specified).
#input: 1. seqFile (headers must contain ">")
#		2. num of chars per seq line (optional)
#output: FASTA file to stdout
#example: perl convertSeqFileFormat.pl input.fa
#example: perl convertSeqFileFormat.pl input.fa 0 0 replace
#example: perl convertSeqFileFormat.pl input.fa embl fasta (didn't try this one yet; should work).
use strict; 
use Bio::SeqIO;

(my $seqFile, my $informat, my $outformat, my $replace) = @ARGV; 

$informat = "fasta" unless $informat; 
$outformat = "fasta" unless $outformat; 

my $inseq = Bio::SeqIO->new(-file   => "<$seqFile",-format => $informat);
my $seq_out = Bio::SeqIO->new(-fh   => \*STDOUT,-format => $outformat);

while (my $seq = $inseq->next_seq) {
    $seq_out->write_seq($seq); 
}








