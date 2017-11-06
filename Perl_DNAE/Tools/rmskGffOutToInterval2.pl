#converts the gff file outputed by RepeatMasker to the Interval file formats
#GFF line example: 
# chr1    RepeatMasker    similarity      1206    1516    17.1    +       .       Target "Motif:L2-1B_ACar" 4724 5049
#INPUT: the gff file (the output file will be same name but .interval instead of .gff)
#Note: the script works for a library where all sequences are from the same Class + Family (all get same annotation)
#Note: appends, so that multiple outputs can be inserted to same interval file (useful when lib contains subfamilies from 
#		different classes/families). 
#version 2 Notes: 	1. Make sure that families don't have slashes: cat RepeatMasker.out.gff | sed 's/\///' |  grep '\/' --- should produce 0 lines
#					2. In general check that parsing is valid for the specific output.
#					3. creates a bed file with the 'name' field which is good for the respective FASTA header
#Modification 7/31/12: The problem: if no Family is listed the original regex failed. This happened for DNA, DNA?, scRNA, and srpRNA. 
#						So I added a if-else to deal with this. 
use strict; 
(my $gff_file, my $assembly) = @ARGV; 
my $class, my $family; 
my $invl_file = $gff_file; 
my $bed_file =  $gff_file; 
$invl_file =~ s/.gff$/.interval/; 
$bed_file =~ s/.gff$/.bed/; 

open(GFF, "<$gff_file") or die "couldn't open $gff_file\n"; 
open(INVL, ">>$invl_file") or die "couldn't open $invl_file\n"; 
open(BED, ">>$bed_file") or die "couldn't open $bed_file\n"; 
while(<GFF>){
	chomp;  
	next if $_ =~ /^\#/; 
	my @fields = split (/\t/, $_); #change: split by tab
	(my $chr, my $start, my $end, my $strand, my $id) = ($fields[0], $fields[3], $fields[4], $fields[6], $fields[8]);
	$start--; #1-base to 0-base start.
	my $subfamily, my $class, my $family; 
	if ($id =~ /Target=(\S+).+Class=(\S+)\/([^; ]+)/){
		($subfamily, $class, $family) = $id =~ /Target=(\S+).+Class=(\S+)\/([^; ]+)/; 
	}
	elsif ($id =~ /Target=(\S+).+Class=([^;\/]+)/){
		($subfamily, $class) = $id =~ /Target=(\S+).+Class=([^;\/]+)/; 
		$family = $class; 
	}
	print  INVL $chr ."\t". $start."\t". $end ."\t". $strand ."\t". $subfamily ."\t". $class ."\t". $family ."\n"; 
	my $name = $assembly."_".$chr."_".$start."_".$end."_".$strand;
	print  BED $chr ."\t". $start."\t". $end ."\t". $name ."\t". "." ."\t". $strand ."\n";
}
close (BED); 
close(GFF);
close(INVL);