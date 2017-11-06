#Function: convert interval file to bed format which can be used to extract sequences from genome with bedtool's getfasta (with -name)
use strict; 
(my $invl_file, my $assembly) = @ARGV; 
my $bed_file = $invl_file; 
$bed_file =~ s/(.interval|.invl)$/.bed/;

open(INVL, "<$invl_file") or die "couldn't open $invl_file\n"; 
open(BED, ">$bed_file") or die "couldn't open $bed_file\n"; 
while(<INVL>){
	chomp;  
	next if $_ =~ /^\#/; 
	my @fields = split (/\s+/, $_); 
	(my $chr, my $start, my $end, my $strand) = ($fields[0], $fields[1], $fields[2], $fields[3]);
	my $name = join ("_",$assembly, $chr, $start, $end, $strand); 
	print BED $chr ."\t". $start."\t". $end ."\t". $name ."\t". "." ."\t". $strand ."\n"; 
}
close(INVL);
close(BED);