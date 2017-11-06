# A package of functions for nucleotide frequencies
###NOTE: this script may be unneeded because of scripts I wrote in analysisSubs.
package nucProbs; 
use strict; 
#use Bio::SeqIO; 


sub freqFromSeq{
	(my $seq, my $nuc, my $flank) = shift; 
	$seq =~ s/\s+//g; #erase whitespaces from seq
	$seq = lc $seq; 
	my @seqAr = split("",$seq);
	my $winSize = $flank * 2 + 1; 
	my @nucs = ("a", "c", "g", "t"); 
	my %map = map {$nucs[$_] => $_} (0..3); #map 'a' to 0, 'c' to 1 etc. 
	print $seq ."\n"; 
	my @freqs = (); 
	foreach my $i (0..($winSize-1)) { foreach my $j (0..3) { $freqs[$i][$j] = 0; } } #lim x flank_size matrix of zeros
	foreach (my $i=$flank; $i<=$#seq-$flank; $i++){
		if ($seqAr[$i] eq $nuc){
			
			
		}
	}
	
	
	my $regex = "([actg]{".$flank."})(".$nuc.")([actg]{".$flank."})"; 
	while ( $seq =~ m/$regex/g ){
		print $& ."\n"; 
	}
	
	
	
	return \@freqs; 
}



sub freqFromFile{

}







return 1; 