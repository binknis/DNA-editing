#1. converts decription line format of sequences' defline downloaded directly from UCSC Genome Browser
#to format retrieved by downloading via Galaxy from UCSC. 
#2. Files from Genome Browser have "1-base" starts, and should be converted to "0-base" for consistency with
#	the files downloaded directly from UCSC.


use strict; 

(my $gb_format, my $galaxy_format) = @ARGV; 
my $start_0_base; 
open (BROWSER, "<$gb_format") || die ("couldn't open $gb_format\n"); 
open (GALAXY, ">$galaxy_format") || die ("couldn't open $gb_format\n"); 
while (my $line = <BROWSER>){
	if ($line =~ /^>/){
		$line =~ /^>([^_]+)_.+range=(\S+):(\d+)-(\d+).+strand=([+-])/; 
		$start_0_base = $3; 
		$start_0_base--; 
		$line = ">" . $1 ."_". $2 ."_". $start_0_base ."_". $4 ."_". $5 ."\n";
	}
	print GALAXY $line; 
}
close(BROWSER); 
close(GALAXY); 