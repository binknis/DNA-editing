use strict; 

#convert format is: 
#from "418=mm9=chr5:138287289-138287516+=LTR=ERVK=IAP-d-int" 
#to "IAPLTR4_I-chr17:19178678-19179524-" format (Shai's). 

(my $clustersFile, my $convertFormat) = @ARGV; 
open (CLUSTERS, $clustersFile) || print "couldn't open $clustersFile!\n"; 
my %editedNames = (); 
my $a_name; 
while (<CLUSTERS>){
	if ($_ =~ /^A name = (\S+)/){
		$a_name = $1; 
		if ($convertFormat){
			my @ar = split(/=/,$a_name); 
			$a_name = $ar[5] ."-". $ar[2]; 
		}
		if (not exists $editedNames{$a_name}){
			$editedNames{$a_name}=0; 
		}
		$editedNames{$a_name}++; 
	}

}
close (CLUSTERS);  
foreach my $key (sort (keys %editedNames)){
	print $key . "\t" . $editedNames{$key} . "\n"; 
}


