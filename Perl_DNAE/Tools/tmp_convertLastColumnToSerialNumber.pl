use strict; 

my $count = 1; 
while (<STDIN>){
	if ($_ =~ /^track/){
		print $_; 
		next; 
	}
	my $line = ($_ =~ /^(.+)\t\d+\n/); 
	print $1 ."\t". $count ."\n"; 
	$count++; 
}