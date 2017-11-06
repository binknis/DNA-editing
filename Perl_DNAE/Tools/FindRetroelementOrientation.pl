
#Runs blastx against GAG, POL and ENV proteins. 
#The result can be checked by the Query coordinates: 
#	if the coordinates increase in the alignment then they're in correct orientation
#	if the coordinates decrease in the alignment then they're in opposite orientation. 

use strict; 
(my $reDNAFile, my $reProtFile, my $outFile) = @ARGV; 

my $blastx = "/usr/common/Software/Blast/blast-2.2.23/bin/blastall -p blastx -d $reDNAFile -i $reProtFile -e 1e-3 -F F > $outFile";
$blastx =~ s/\(/\\\(/g; $blastx =~ s/\)/\\\)/g; #add backslashes before "(" and ")" symbols. 
system($blastx);
