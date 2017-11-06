#creates a table from tuples, using 3 fields in each tuple: x axis, y axis, and "value". 
use strict; 
(my $tupleFile, my $xInd, my $yInd, my $valInd) = @ARGV; 
my %pval_th_pairs = (); 
open (TUP, $tupleFile) || die "couldn't open clustStatFile\n"; 
while (my $line = <TUP>){
	chomp; 
	my @fields = split (/\s+/, $line); 
	my $val = $fields[$valInd]; 
	$val =~ s/%//; 
	$pval_th_pairs{$fields[$xInd]}{$fields[$yInd]} = $val; 
}
close(TUP); 

foreach my $pval (sort{$a <=> $b}(keys (%pval_th_pairs))){
	foreach my $th (sort{$a <=> $b}(keys (%{$pval_th_pairs{$pval}}))){
		print $pval_th_pairs{$pval}{$th} . " ";
	}
	print "\n"; 
}