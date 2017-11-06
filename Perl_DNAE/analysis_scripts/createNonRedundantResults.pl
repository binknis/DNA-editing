# Function: removes redundancy from track output
# Input: Data needed to access a specific Tracks directory
# Output: 
# Notes: 1. calls perl516

use File::Path qw(mkpath);
use strict;


( my $organism, my $class, my $family, my $pval, my $th, my $control, my $subdir ) = @ARGV;
$pval = "1e-" . $pval if $pval =~ /^\d+$/;
$subdir = 0 if $subdir eq ''; 
$control = 0 if $control eq '';
my $pmotif = '1e-3'; 

my $dir = "../Data/" . $organism . "/" . $class . "/results"; 
$dir .= "/$subdir" unless ($subdir eq 0 or $subdir eq ""); 
$dir .= "/Tracks/tracks"; 
my $suffix = "_" . $organism . "_" . $class . "_" . $family . "_" . $pval . "_" . $th . ($control ? "_control" : ""); 
$dir .= $suffix; 

#create track files if doesn't exist (avoids the need to run two commands
if (not -d $dir){
	system("perl516 analysis_scripts/createTrackFiles.pl $organism $class $family $pval $th $control $subdir"); 
}




