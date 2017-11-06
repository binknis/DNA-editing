#used to convert the fasta output of fastaToBed to the input fasta format needed for sortGenome.pl
#changes: 1. adds assembly name at head of each defline. 2. replaces spacers with underscores.

use strict; 
use File::Copy;
(my $fasta, my $assembly) = @ARGV; 


open (FASTA, "<$fasta") || die "couldnt open $fasta\n"; 
while (my $line = <FASTA>){
	chomp $line; 
	my $genome; 
	my $pos;
	if ($line =~ /^>/){
		$line =~ s/\)//; $line =~ s/\(//; 
		if ($line =~ /^>([^:]+):(\d+)-(\d+)([+-])\s*$/){ #original fasta format from UCSC output
			$line = ">".$assembly."_".$1."_".$2."_".$3."_".$4; 
		}
		elsif ($line =~ /^>/){
			die "bad defline format\n"; 
		}
	}
	print $line ."\n";
}
close (FASTA); 