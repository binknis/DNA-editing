#Function: Finds existence and integrity of ORFs a DNA sequence (typically retroelement/virus)
#Input: 1. protFile: Fasta file of proteins to search for in the sequence
use strict; 
use Bio::SeqIO;

my $seqio  = Bio::SeqIO->new( '-format' => 'embl' , -file => 'myfile.dat');
while ($seqobj = $seqio->next_seq()){
	my $protSeq = $seqobj->translate; 
}