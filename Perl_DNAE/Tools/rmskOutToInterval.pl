#converts the rsmk.out file outputed by RepeatMasker to the Interval file formats
#rsmkout line example: 
#1964   13.8  3.4  0.6  Scaffold14911    6180   6506   (1982) C lcl|CR1-Y_Pass        LINE/CR1               (6)   3782    3447    1
#INPUT: the rmsk.out file (the output file will be same name but .interval instead of .RepeatMasker.out)

#Notes: 	1. Make sure that families don't have slashes: cat RepeatMasker.out.gff | sed 's/\///' |  grep '\/' --- should produce 0 lines
#			2. In general check that parsing is valid for the specific output.
#			3. creates a bed file with the 'name' field which is good for the respective FASTA header

use strict; 
(my $rmsk_file, my $assembly, my $dropAsteriskMarked, my $invlOutdir, my $bedOutdir) = @ARGV; 
my $invl_file = $rmsk_file; 
my $bed_file =  $rmsk_file; 
$invl_file =~ s/.RepeatMasker.out$/.interval/; 
$bed_file =~ s/.RepeatMasker.out$/.bed/; 
if($invlOutdir){ #set output dir (if not default = same as rmsk)
	$invl_file =~ s/\S+\//$invlOutdir\//;
}
if($bedOutdir){ #set output dir (if not default = same as rmsk)
	$bed_file =~ s/\S+\//$bedOutdir\//;
}

open(RMSK, "<$rmsk_file") or die "couldn't open $rmsk_file\n";
open(INVL, ">$invl_file") or die "couldn't open $invl_file\n";
open(BED, ">$bed_file") or die "couldn't open $bed_file\n";
while(my $l = <RMSK>){
	chomp $l;  
	next if $l =~ /^(\#|\s*(SW|score))|^\s*$/; #skip comments, headers and empty lines
	if($dropAsteriskMarked){
		next if $l =~ /\*\s*$/; 
	}
	$l =~ s/\s+/\t/g;
	$l =~ s/^\s+//;
	
	my @fields = split (/\t/, $l); #change: split by tab
	(my $chr, my $start, my $end, my $strand, my $subfam_raw, my $classFam) = ($fields[4], $fields[5], $fields[6], $fields[8], $fields[9], $fields[10]);
	
	#modify input
	$start--; #1-base to 0-base start
	$strand =~ s/C/-/; #strand format change
	(my $subfam) = $subfam_raw =~ /([^|]+)$/; #extract 'name' (subfam)
	(my $class, my $fam) = split('/', $classFam); 
	
	#print interval and BED files
	print  INVL $chr ."\t". $start."\t". $end ."\t". $strand ."\t". $subfam ."\t". $class ."\t". $fam ."\n"; 
	my $name = $assembly."_".$chr."_".$start."_".$end."_".$strand;
	print  BED $chr ."\t". $start."\t". $end ."\t". $name ."\t". "." ."\t". $strand ."\n";
}
close (BED); 
close(RMSK);
close(INVL);
