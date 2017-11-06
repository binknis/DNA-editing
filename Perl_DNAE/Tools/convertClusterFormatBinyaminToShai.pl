#FUNCTION: converts a cluster file in Binyamin's format to Shai's format. 
#			used to check that algorithm produces similar results. 
use strict; 

# Found cluster:
# G name = 74=mm9=chr1:145680148-145681649-=LTR=ERVK=IAPEY3-int
# A name = 889=mm9=chr9:89939645-89941263-=LTR=ERVK=IAPEY3-int
# Total probability = 5.64040476510701e-20, 1 cluster
# Edited MM serial no.    :     8     9    10    11    12    13    14    15    16    17    18    19    20    21    22
# Locations G seq:  1260  1266  1284  1287  1293  1296  1307  1313  1315  1325  1368  1382  1387  1405  1410
# Locations A seq:  1379  1385  1403  1406  1412  1415  1426  1432  1434  1444  1487  1501  1506  1524  1529
# Total length: 421, Direct mismatches: 24, All mismatches: 29
# End cluster

# Found cluster:
# A name = IAPLTR4_I-chr17:19178678-19179524-
# G name = IAPLTR4_I-chr2:11865258-11865620-
# Total probability = 8.04092854440382e-13, 1 cluster
# Edited MM serial no.    :    11    12    13    14    15    16    17    18    19    20    21    22    23
# All locations (in A seq):   205   208   209   232   240   244   250   258   262   269   274   289   291
# Total length: 311, Direct mismatches: 29, All mismatches: 39
# End cluster


my $clustFile = shift; 

open (SHAI, "<$clustFile") || die "$clustFile didn't open\n"; 
open (BK, ">".$clustFile."_convertedToBinyaminFormat") || die "New cluster file didn't open\n"; 

while (my $line = <SHAI>){
	#reached line that is different in the formats
	if ($line =~ /^G name = (\S+)/){
		#get G label line
		my $G_label_line = convertLabelLine($line); 
		#get A label line
		$line = <SHAI>; 
		my $A_label_line = convertLabelLine($line);
		#print in reverse order
		print BK $A_label_line; 
		print BK $G_label_line; 
		#print two lines that don't need to be changed
		print BK ($line = <SHAI>); 
		print BK ($line = <SHAI>); 
		#skip G seq positions (not in Shai's files).
		$line = <SHAI>; 
		#modify and print A positions
		$line = <SHAI>; 
		(my $locations) = ($line =~ /^Locations A seq:(.+)$/); 
		print BK "All locations (in A seq):" . $locations."\n";
	}
	else{
		print BK $line; 
	}
}

close (SHAI); 
close (BK); 

sub convertLabelLine{
	my $old = shift; 
	(my $new, my $oldLabel) = ($old =~ /^([AG] name = )(\S+)/); 
	my @fields = split('=',$oldLabel); 
	$new .= $fields[5] ."-". $fields[2] . "\n"; 
	return $new; 
}