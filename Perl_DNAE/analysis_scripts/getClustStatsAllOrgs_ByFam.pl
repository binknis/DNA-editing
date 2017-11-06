use strict; 
use getAll; 

(my $class, my $fam) = @ARGV; 

my $orgs = getAll::organisms(); 
foreach my $org (@$orgs){
	my $clust_stat_file = getAll::clustStatsFamAsFilename($org, $class, $fam, 1);
	next unless $clust_stat_file;
	my $lines = getAll::lines($clust_stat_file);
	
	foreach my $line (@$lines){
		chomp $line;
		print "$org\tLINE\t$line\n";
	}
}