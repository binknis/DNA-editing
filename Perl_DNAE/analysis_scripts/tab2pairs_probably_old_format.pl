#Function: converts tab delimited pairs to coords format

use strict; 

while (my $l = <>){
	chomp $l; 
	my @f = split (/\t/, $l); 
	print $f[0] .":". $f[1] ."-". $f[2] . $f[3] ."\t". $f[4] .":". $f[5] ."-". $f[6] . $f[7] ."\n"; 
}

