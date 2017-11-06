use strict; 


(my $infile, my $outfile) = @ARGV; 

open(FH, $infile) or die "open $infile failed!";
open(OUT, ">".$outfile) or die "open $outfile failed!";
my $is_seq = 0; 
while(my $l = <FH>){
	chomp $l;
	if($l =~ /^@/){
		$l =~ s/@(\S+).*/>\1/; 
		print OUT $l ."\n"; 
		$is_seq=1; 
		next; 
	}
	elsif($l =~ /^\+/){
		$is_seq=0; 
	}
	
	if($is_seq==1){
		print OUT $l ."\n"; 
	}
}
close(FH); 
close(OUT); 



