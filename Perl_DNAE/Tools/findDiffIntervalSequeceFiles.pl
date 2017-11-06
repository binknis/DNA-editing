#INPUT: interval file, sequence file and Organism name. 
#OUTPUT: 	1.all sequences that are present in interval file and not in sequence file. 
#			2. Visa versa. 
#The output is outputted to STDOUT
use strict; 
(my $seqDataFile, my $seqFile, my $organism) = (@ARGV);
my %coords = (); 
my @extraIntervals = (); 
my @extraSeqs = (); 

#get coordinates from sequence and interval files
open(INTERVAL,"<$seqDataFile") || die ("couldn't open $seqDataFile file\n"); 
while (my $line = <INTERVAL>){
	my @fields = split(/\s+/, $line); 
	my $oneCoord = $fields[0]."_".$fields[1]."_".$fields[2]."_".$fields[3]; 
	if (exists $coords{$oneCoord}) {
		$coords{$oneCoord}++;
	}
	else{
		$coords{$oneCoord}=1; 
	}
}
close(INTERVAL); 

open(SEQ,"<$seqFile") || die ("couldn't open $seqFile file\n"); 
while (my $line = <SEQ>){
	if ($line =~ /^>/){
		$line =~ /[^_]+_(\S+)/;
		if (exists $coords{$1}) {
			$coords{$1} -= 1;
		}
		else{
			$coords{$1} = -1; 
		}
	}
}
close(SEQ); 

#check for redundancy
foreach my $oneCoord(keys(%coords)){
	if ($coords{$oneCoord} > 0){
		push(@extraIntervals,$oneCoord);
	}
	elsif($coords{$oneCoord} < 0){
		push(@extraSeqs,$oneCoord);
	}
}

print "Total extra Intervals:\t" . $#extraIntervals+1 . "\n"; 
print "Total extra Sequences:\t" . $#extraSeqs+1 . "\n"; 
print "Extra intervals list: \n"; 
foreach my $oneCoord(@extraIntervals){
	print $oneCoord . "\n"; 
}
print "Extra sequences list: \n"; 
foreach my $oneCoord(@extraSeqs){
	print $oneCoord . "\n"; 
}