#converts the gff file outputed by RepeatMasker to the Interval file formats
#GFF line example: 
# chr1    RepeatMasker    similarity      1206    1516    17.1    +       .       Target "Motif:L2-1B_ACar" 4724 5049
#INPUT: the gff file (the output file will be same name but .interval instead of .gff)
#Note: the script works for a library where all sequences are from the same Class + Family (all get same annotation)
#Note: appends, so that multiple outputs can be inserted to same interval file (useful when lib contains subfamilies from 
#		different classes/families). 
use strict; 
(my $gff_file, my $class, my $family) = @ARGV; 
my $invl_file = $gff_file; 
$invl_file =~ s/.gff$/.interval/; 

open(GFF, "<$gff_file") or die "couldn't open $gff_file\n"; 
open(INVL, ">>$invl_file") or die "couldn't open $invl_file\n"; 
while(<GFF>){
	chomp;  
	next if $_ =~ /^\#/; 
	my @fields = split (/\s+/, $_); 
	(my $chr, my $start, my $end, my $strand, my $name) 
						= ($fields[0], $fields[3], $fields[4], $fields[6], $fields[9]);
	$start--; #1-base to 0-base start.
	$name =~ s/(\"|Motif:)//g; 
	print  INVL $chr ."\t". $start."\t". $end ."\t". $strand ."\t". $name ."\t". $class ."\t". $family ."\n"; 
}
close(GFF);
close(INVL);
