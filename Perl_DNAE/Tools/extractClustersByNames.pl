#This script must be run from one level down from "Data"
#creates a new folder inside results of Data/Organism/Class and creates new clusters inside it.


#Reminder: label format: 1089=hg19=chr1:145376886-145377097-=LINE=L1=HAL1

use strict; 
use File::Path qw(mkpath); 
use getAll; 
# perl Tools/extractClustersByNames.pl Mouse LTR ERVK 0 MusD RLTRETN_Mm ETnERV2-int ETnERV3-int ETnERV-int MMETn-int
(my $organism, my $class, my $family, my $subdir, my $outName, my @names) = @ARGV;

#insert names into hash
my %names = (); 
foreach my $argName (@names){
	$names{$argName}=1; 
}

#construct directory paths (raw and new clusters)
my $cluster_dir = "../Data/$organism/$class/results"; 
$cluster_dir .= "/$subdir" if $subdir != 0;
my $new_cluster_dir = $cluster_dir ."/Clusters_". $outName; 

#get cluster file names
my $clusterFiles = getAll::clustersFam($organism, $class, $family, 0, $subdir); 
die "no cluster files obtained for $outName\n" if ($clusterFiles == 0); #no clusters obtained
mkpath($new_cluster_dir); 

# (file names will be same as in old directory, except for the family which will be changed to $outName) 
#copy and extract from each file into new file
foreach my $clustFile (@$clusterFiles){
	open (OLD, "<$cluster_dir/$clustFile") || die "$cluster_dir/$clustFile didn't open!\n";
	$clustFile =~ /clusters_([^_]+)_([^_]+)_(\S+)_(1e-\d+)_(\d+)(_control)?.txt$/; 
	my $clustOutFile  = "clusters_".$1."_".$2."_".$outName."_".$4."_".$5.$6.".txt";
	open (NEW, ">$new_cluster_dir/$clustOutFile") || die "$new_cluster_dir/$clustOutFile didn't open!\n";
	while (my $line = <OLD>){
		if ($line =~ /^Found cluster:/){
			my $notInNameList = 0; 
			my $cluster = $line; 
			while ($line = <OLD>){
				if ($line =~ /^[GA] name = (\S+)/){
					(my $name) = ($line =~ /^[GA] name = \S+=(\S+)/); #get name
					$notInNameList = 1 unless exists $names{$name}; 	
				}
				$cluster .= $line;
				last if ($line =~ /^End cluster/); 
			}
			print NEW $cluster unless ($notInNameList);
		}
	}
	close(NEW); 
	close(OLD);
}