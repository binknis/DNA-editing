#prints frequencies of lengths of sequences in a cluster file
#counts A and/or G sequences. 
#INPUT: 1. cluster file - full path; 2. G, A or GA - as a flag of which sequences' lengths to fetch. 
#OUTPUT: Length and frequency - to STDOUT. 

use strict; 

(my $clusterFile, my $GA_flag) = @ARGV; 
my $ga; 
my %lens; 
if ($GA_flag == 1){
	$ga = "G"; 
}
elsif($GA_flag == 2){
	$ga = "A";
}
else {
	$ga = "GA";
}


open (CLUST, $clusterFile) || die "didn't open $clusterFile\n"; 
while (my $line = <CLUST>){
	if ($line =~ /[$ga] name/){
		#print $line; #***
		chomp $line; 
		($line) = $line =~ /[GA] name = (\S+)/;
		my @fields = split (/=/,$line); 
		(my $start, my $end) = $fields[2] =~ /:(\d+)-(\d+)/;
		my $len = $end - $start; 
		$lens{$len} = 0 unless exists $lens{$len}; 
		$lens{$len}++;
	}
}
close (CLUST); 

foreach my $length (sort keys %lens){
	print $length . "\t" . $lens{$length} . "\n"; 
}