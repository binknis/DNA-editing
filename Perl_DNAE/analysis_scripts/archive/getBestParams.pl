use strict; 

use lib 'analysis_scripts'; 
use getParams; 

(my $org, my $class, my $family, my $subdir, my $target_val, my $margin) = @ARGV; 
my $cstats_file = "../Data/$org/$class/results/$subdir/cluster_stats_".$org."_".$class."_".$family.".txt"; 

my $closest_to_col = 8; # 8 - FP edited
my $maximize_by_col = 6; # 6 - # edited (real; G-to-A)

print $org ."\t". $class ."\t". getParams::closestToMaximizeBy_returnLine($cstats_file, $closest_to_col, $target_val, $maximize_by_col, $margin);


