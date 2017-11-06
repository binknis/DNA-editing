#Function: renames all cluster_stats_ files so they include the Organism and Class (in addition to the family). 
use strict; 
use getAll; 

my $organisms = getAll::organisms(); 
my $resDir; 
my @resFiles = (); 
foreach my $org (@$organisms){
	my $classes = getAll::classes($org); 
	foreach my $class (@$classes){
		$resDir = "../Data/$org/$class/results";
		@resFiles = glob $resDir . "/cluster_stats_*";
		foreach my $oldName (@resFiles){
			my $newName = $oldName; 
			my $newFilePrefix = "cluster_stats_".$org."_".$class."_"; 
			$newName =~ s/cluster_stats_/$newFilePrefix/; 
			rename $oldName, $newName; 
		}
	}
}