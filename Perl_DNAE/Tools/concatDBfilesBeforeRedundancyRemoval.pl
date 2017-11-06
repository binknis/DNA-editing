#Function: concatanates the content of DB files from one org into the DB files in another org. 
#			The next step is to remove redundancy of the target DB files, with removeRedundancyFromDBfiles.pl
use strict; 
use getAll;
use File::Copy;  
(my $org_from, my $org_to, my $class) = @ARGV; 
my @prefixes = ("LenStats", "Len", "Nuc");
my $fams = getAll::families($org_from, $class); 
for my $pref(@prefixes){
	for my $fam(@$fams){
		my $from_file = "../Data/$org_from/$class/db/".$pref ."_". $fam .".txt";
		my $to_file = "../Data/$org_to/$class/db/".$pref ."_". $fam .".txt";
		system ("cat $from_file >> $to_file");
	}
}