#Function calculates diff in output between db- and blast-count files created by validateCompletionOfAllBlasts.pl
#Note: I used this by concatenating all output files from all validateCompletionOfAllBlasts.pl runs and it outputed for all organisms-classes-families
use strict; 
use lib $ENV{HOME} ."/Perl_DNAE"; 
use getAll; 
(my $db_file, my $blast_file) = @ARGV; 
my %db = (); 
my %diff = (); 
my %blast = (); 
open (DB, $db_file) || die "open $db_file\n"; 
while (my $line = <DB>){
	chomp $line; 
	(my $taxa, my $num) = $line =~ /^(.*)\s+(\d+)/; 
	$db{$taxa} = $num; 
	$diff{$taxa} = $num; 
}
close (DB);
 
open (BLAST, $blast_file) || die "open $blast_file\n"; 
while (my $line = <BLAST>){
	chomp $line; 
	(my $taxa, my $num) = $line =~ /^(.*)\s+(\d+)/; 
	$blast{$taxa} = $num; 
	$diff{$taxa} = $db{$taxa} - $num; 
}
close (BLAST);

print "DB count minus BLAST count:\n"; 
foreach my $taxa (sort keys %db){
	print $taxa ."\t". $db{$taxa} ." - ". $blast{$taxa} . " = ". $diff{$taxa} . "\n" if $diff{$taxa} != 0; 
}