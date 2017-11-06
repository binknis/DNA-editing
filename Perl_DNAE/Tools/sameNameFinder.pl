#Function: "Same name, different Class or family" finder. 
#INPUT: interval file(downloaded from Galaxy, as in Raw_data directory)
use strict; 
( my $seqDataFile) = @ARGV;

### Create "familyTree", i.e. tree describing the hierarchy of Class/Family/Name of all sequences ###
open(SEQDATA,"<$seqDataFile") || die ("couldn't open description file"); 
my %nameToTaxa = (); #chr_start_end => Class=Family=Name
my $class; 
my $family; 
my $name; 

while (my $line = <SEQDATA>){
	chomp; 
	next if ($line =~ /^\s*#/); #skip remark lines
	
	my @seqData = split(/\s+/,$line); 
	$class = $seqData[5];
	$family = $seqData[6]; 
	$name = $seqData[4]; 

	#replace slashes with underscores; remove question-marks. (format used in sortGenome.pl).
	foreach my $taxa($class, $family, $name){
		$taxa =~ s/\//_/g; 
		$taxa =~ s/\?$//;
	}
	
	$nameToTaxa{$name}{"$class=$family"}=0; 
}
close(SEQDATA); 

foreach $name (keys(%nameToTaxa)){
	if (scalar(keys(%{$nameToTaxa{$name}})) > 1){
		print "$name:\n"; 
		foreach my $taxa(keys(%{$nameToTaxa{$name}})){
			print "$taxa\t"; 
		}
		print "\n\n"; 
	}
}