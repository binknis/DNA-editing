#Function: convert bed or gff file to ANNOVAR input
#*** under construction. 
use strict; 


(my $from_to, my $inFile, my $outFile) = @ARGV; 
(my $ref, my $mut) = split('', uc $from_to); 

my $in; 
my $out; 
if (not $inFile || $inFile =~ /STDIN/){ #write to stdin
	open ($in, \*STDIN) or die "STDIN didn't open\n"; 

}
else{

}

while (my $l = <$in>){

}

if ($outFile){

}
else{

}


close ($out) if ($outFile); 

