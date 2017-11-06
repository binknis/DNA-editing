#B must fully contain any A sequence and prints B.
#example: A are exapted elemnts, B are edited elements. prints all edited containing an expated element
use strict; 
(my $a, my $b) = @ARGV; 

my %a_hash = (); 
open (A, "<".$a) or die "$a didn't open\n"; 
while (my $a_coords = <A>){
	(my $chr, my $start, my $end) = split (/\t/, $a_coords);
	$a_hash{$chr} = () unless exists $a_hash{$chr}; 
	push (@{$a_hash{$chr}},$start ."-". $end);
}
close (A); 

open (B, "<".$b) or die "$a didn't open\n"; 
while (my $b_coords = <B>){
	chomp $b_coords; 
	(my $chr_b, my $start_b, my $end_b) = split (/\t/, $b_coords); 
	foreach my $a_coords(@{$a_hash{$chr_b}}){
		(my $start_a, my $end_a) = split (/-/, $a_coords); 
		if ($start_a >= $start_b && $end_a <= $end_b){
			# print $b_coords ."\t".$start_a."\t".$end_a."\n";
			print $b_coords ."\n";
			last;
		}
	}
}
close(B); 
