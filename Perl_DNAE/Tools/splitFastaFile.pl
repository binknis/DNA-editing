#split a fasta file into multiple files with (approximately) the same amount of sequences
use strict;
use Bio::SeqIO;


(my $in_file, my $out_prefix, my $num_files) = @ARGV;

#get num seqs in input file
my $seqNum = `grep -c ">" $in_file`;
my $seqPerFile = int($seqNum / $num_files); 
if ($seqPerFile == 0){ #less sequences than output files demanded - use only one output file
	$num_files = 1; 
	$seqPerFile = $seqNum; 
}

#create output streams
my %out = (); 
for (my $i=1; $i<=$num_files; $i++){
	$out{$i} = new Bio::SeqIO(-file => ">$out_prefix.$i", -format=>'fasta');
}

#read input and write output
my $in  = new Bio::SeqIO(-file  => $in_file);
my $i=1; 
my $seqOutCount = 0; 
while (my $seq = $in->next_seq) {
	$out{$i}->write_seq($seq);
	$seqOutCount++;
	unless ($seqOutCount < $seqPerFile || $i == $num_files){ #if enough seqs per file were written - start writing in next file (unless this is last file)
		$i++; 
		$seqOutCount = 0;
	}
}

