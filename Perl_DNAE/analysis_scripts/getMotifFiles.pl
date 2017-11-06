use strict; 
my @orgs = qw/Baboon	Bushbaby	Chimp	Gibbon	Gorilla	Human	Marmoset	MLemur	Orangutan	Rhesus	SMonkey	Tarsier/; 
foreach my $org (@orgs){
	print $org ."_". "ERV1_1e-7_7"."\n"; 
	system ("cat $org/LTR/results/Tracks/tracks_".$org."_LTR_ERV1_1e-7_7/logo_G_".$org."_LTR_ERV1_1e-7_7_freq.txt")
}
