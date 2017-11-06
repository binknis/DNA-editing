use strict; 
use getAll; 

(my $class, my $family, my $outPath) = @ARGV; 

my $cluster_stats_files = getAll::clustStatsFamAllOrgs($class, $family, 1); 
open (OUT, ">$outPath") || die "didnt open $outPath\n"; 
foreach my $csf (@$cluster_stats_files){
	(my $org) = $csf =~ /\/Data\/([^\/]+)\//; #get organism name from path
	# print $org ."\n";
	my $lines_ref = getAll::lines($csf); 
	foreach my $line (@$lines_ref){
		chomp $line; 
		print OUT $org ."\t". $class ."\t". $line ."\n";
	}
}

