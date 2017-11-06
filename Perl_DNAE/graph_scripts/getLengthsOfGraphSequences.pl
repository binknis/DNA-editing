
#output flag: 0: length frequency. 1: min and max length

use strict; 
use Bio::SeqIO;
use List::Util qw[min max];

(my $graphFile, my $outputFlag) = @ARGV; 
my %nameToSerials = (); 
my %lenHist = ();
my $serial; 
my $len; 

#open graph file and create a hash-of-hash of Elements => Name => Sequence_serial_number
open (GRAPH, "<$graphFile") || die "couldn't open $graphFile\n"; 
while (my $line = <GRAPH>){
	chomp $line; 
	my @sequences = split(/\s+/,$line);
	$sequences[0] =~ /^([^=]+)=.+=([^=]+)$/; #source sequence
	$nameToSerials{$2}{$1} = 1; 	#{$name}{$serial}=1
	$sequences[1] =~ /^([^=]+)=.+=([^=]+)$/;
	$nameToSerials{$2}{$1} = 1; #target sequence
}
close(GRAPH); 
#open each of the Name files and count each length
$graphFile =~ /graph_([^_]+)_([^_]+)_([^_]+)_(1e-\d+)_(\d+)(_control)?\.txt$/;
my $seqDir = "../Data/$1/$2/db/files_".$3; 
unless (-d $seqDir) {print $seqDir . " doesn't exist!\n"; exit;}
foreach my $name (sort(keys(%nameToSerials))){
	my $nameFile = $seqDir . "/Seq_" . $name; 
	unless (-f $nameFile){ #name file doesn't exist
		print $nameFile . " doesn't exits!\n"; 
		next; 
	}
	#read each sequence and insert their length into a hash which holds a histogram of the lengths
	my $seqio_obj = Bio::SeqIO->new(-file => $nameFile);
	while (my $seq_obj = $seqio_obj->next_seq){
		$seq_obj->display_id =~ /^(\d+)=/; 
		$serial = $1;
		next unless exists $nameToSerials{$name}{$serial}; #skip sequences that weren't in the input file
		$len = $seq_obj->length;
		$lenHist{$name}{$len} = 0 unless exists $lenHist{$name}{$len};
		$lenHist{$seq_obj->length}++;
	}
}

#print output
foreach my $name (sort(keys(%nameToSerials))){
	if ($outputFlag){
		print "$name: Min: " . min(keys(%{$lenHist{$name}})) . " Max: " . max(keys(%{$lenHist{$name}})) . "\n"; 
	}
	else{
		foreach my $length (sort{$a <=> $b}(keys(%{$lenHist{$name}}))){
			print $length . "\t" . $lenHist{$name}{$length} . "\n"; 
		}
	}
}

