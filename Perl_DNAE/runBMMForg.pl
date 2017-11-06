use strict;
use getAll; 
use lib "/home/alu/binknis/Perl_DNAE/analysis_scripts"; 
use analysisSubs;

(my $org, my $class, my $out) = @ARGV; 

analysisSubs::calcBMMForg($org, $class, $out); 
