use strict; 

use lib 'analysis_scripts'; 
use getParams; 

(my $org, my $class, my $family, my $subdir, my $target_val, my $margin) = @ARGV; 
my $cstats_file = "../Data/$org/$class/results" . ( ($subdir) ? "/$subdir" : "" ) . "/cluster_stats_".$org."_".$class."_".$family.".txt"; 

my $closest_to_col = 10; # 10 - FP edited
my $maximize_by_col = 8; # 8 - # edited (real; G-to-A)

my $line = getParams::closestToMaximizeBy_returnLine($cstats_file, $closest_to_col, $target_val, $maximize_by_col, $margin);
chomp $line; 
if ($line =~ /NA$/){
	print $org ."\t". $class ."\t". $family ."\t". "NA" . "\n";
}
else {
	print $line ."\n"; 
}


# my $line = getParams::inMarginMaximizeBy_returnLine($cstats_file, $closest_to_col, $target_val, $maximize_by_col, $margin);
# chomp $line; 
# if ($line =~ /NA$/){
	# print $org ."\t". $class ."\t". $family ."\t". "NA" . "\n";
# }
# else {
	# print $line ."\n"; 
# }