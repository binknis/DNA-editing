#function: checks if a org and class has a cluster output file with a specific parameter pair 
#		Used to check the ratio of params runned for the org's class
#To run for all orgs: perl foreachOrgClassExec.pl "tmp/findClassesWithoutSpecificParameter.pl xxx yyy 1e-16 16"
#Input: $pval with or w/o '1e-'
#Note: works only for 5 classes 
#looks for the file also in the Clusters_Raw subdir in 'results' dir. 
use strict; 
use getAll; 
(my $org, my $class, my $pval, my $th, my $perClassOrFam) = @ARGV;
exit unless $class =~ /(LINE|LTR|SINE|RC|DNA)/; 
my $fams = getAll::families($org, $class); 
die "No fams for: $org $class\n" if $fams eq 0; 

# my @pvals = ($pval_l .. $pval_h); 
# my @ths = ($th_l .. $th_h); 

$pval =~ s/1e-//; 

if ($perClassOrFam eq "class"){
	my $exists = 1; 
	foreach my $fam (@$fams){
		my $cluster_file = "../Data/".$org."/".$class."/results/clusters_".$org."_".$class."_".$fam."_1e-".$pval."_".$th.".txt"; 
		my $cluster_file_in_raw_dir = "../Data/".$org."/".$class."/results/Clusters_Raw/clusters_".$org."_".$class."_".$fam."_1e-".$pval."_".$th.".txt"; 
		$exists = 0 unless (-e $cluster_file || -e $cluster_file_in_raw_dir);
	}
	print $org ."\t". $class ."\t". $pval. "\t". $th ."\n" unless $exists; 
}
elsif ($perClassOrFam eq "family"){
	foreach my $fam (@$fams){
		my $cluster_file = "../Data/".$org."/".$class."/results/clusters_".$org."_".$class."_".$fam."_1e-".$pval."_".$th.".txt"; 
		my $cluster_file_in_raw_dir = "../Data/".$org."/".$class."/results/Clusters_Raw/clusters_".$org."_".$class."_".$fam."_1e-".$pval."_".$th.".txt"; 
		unless (-e $cluster_file || -e $cluster_file_in_raw_dir){
			print $org ."\t". $class ."\t". $fam ."\t". $pval. "\t". $th ."\n"; 
		}
	}
}
else{
	print "Last parameter must be 'class' or 'family'\n"; 
}