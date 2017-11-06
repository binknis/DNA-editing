#converts a GFF file to BED format
#GFF line example: 
#INPUT: the gff file (the output file will be same name but .bed instead of .gff)
#Note: appends to output file 

use strict; 
(my $gff_file, my $getGene) = @ARGV; 
my $bed_file = $gff_file; 
$bed_file =~ s/.gff$/.bed/; 

open(GFF, "<$gff_file") or die "couldn't open $gff_file\n"; 
open(BED, ">>$bed_file") or die "couldn't open $bed_file\n"; 
while(<GFF>){
	chomp;
	next if $_ =~ /^(\#|track)/ ;
	my @fields = split (/\s+/, $_); 
	(my $chr, my $start, my $end, my $strand) = ($fields[0], $fields[3], $fields[4], $fields[6]);
	$start--; #1-base to 0-base start
	my $name = ".";
	if($getGene){
		($name) = $fields[8] =~ /;gene=([^;]+);/; #get gene name
	}
	else{
		$name = $chr.":".$start."-".$end.$strand; 
	}
	print BED $chr ."\t". $start."\t". $end ."\t". $name ."\t". "." ."\t". $strand ."\n"; 
}
close(GFF);
close(BED);
