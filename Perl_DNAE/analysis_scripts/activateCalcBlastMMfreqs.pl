use strict; 
use FindBin;
use lib "$FindBin::Bin/..";  # use the parent directory of analysis_scripts directory from where this is run
use lib "$FindBin::Bin"; #because this file is in analysis_scripts, this includes that directory
use getAll; 
use analysisSubs; 

(my $blastFile, my $outPerHSPfile, my $outPerMMfile, my $coordsOut, my $append, my $addStartEndStrandFromHSP) = @ARGV; 
analysisSubs::calcBlastMMfreqs($blastFile, $outPerHSPfile, $outPerMMfile, $coordsOut, $append, $addStartEndStrandFromHSP); 