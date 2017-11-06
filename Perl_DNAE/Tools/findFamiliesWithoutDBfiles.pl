#Function: gets a organism and class name and checks that each family 
#			has all the db files. (used to check that the new 'sortGenome.pl' was 
#			used for all organisms and classes. 
#script to use to check all organisms:  
#	perl foreachOrgClassExec.pl "perl Tools/findFamiliesWithoutDBfiles.pl xxx yyy"
use strict; 
use getAll; 

(my $org, my $class) = @ARGV; 
my $fams = getAll::families($org, $class); 
die "No fams for: $org $class\n" if $fams eq 0; 

my $dbdir = "../Data/$org/$class/db";
foreach my $fam (@$fams){
	my $len = $dbdir ."/Len_".$fam.".txt"; 
	my $lenStats = $dbdir ."/LenStats_".$fam.".txt"; 
	my $nuc = $dbdir ."/Nuc_".$fam.".txt"; 
	print "missing: $org $class $fam\n" unless ( -e $len and -e $lenStats and -e $nuc); 
}