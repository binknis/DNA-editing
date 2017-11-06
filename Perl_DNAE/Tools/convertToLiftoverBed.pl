#Function: converts bed/gff/interval file to format good for liftover (will contain original coordinates in output)
#Input: 1. bedFile: bed/gff/interval file (suffix must be exactly one of those)
#		2. liftFile (optional): output file. default - suffix is replaced with "forLift.bed" suffix.
use strict; 
(my $bedFile, my $liftFile) = @ARGV; 

#get type of input file

die "bad input file suffix - must be bed, gff or intreval\n" unless $bedFile =~ /(bed|gff|interval)$/; 
(my $f_type) = $bedFile =~ /(bed|gff|interval)$/; 

unless ($liftFile){
	$liftFile = $bedFile; 
	$liftFile =~ s/$f_type$/forLift.bed/; 
}

open (BED, $bedFile) or die "open $bedFile\n";
open (LIFT, ">" .$liftFile) or die "open $liftFile\n";
while (my $l = <BED>){
	next if $l =~ /^(track|#)/; #skip description lines
	chomp $l;
	my @f = split(/\t/, $l);
	if ($f_type eq 'bed'){
		my $desc = $f[0] .":".$f[1] ."-".$f[2].$f[5];
		print LIFT $f[0] ."\t". $f[1] ."\t". $f[2] ."\t".$desc."\t"."."."\t".$f[5]."\n";
	}
	elsif ($f_type eq 'gff'){
		my $start = $f[3] - 1; 
		my $desc = $f[0] .":". $start ."-".$f[4].$f[6];
		print LIFT $f[0] ."\t". $start ."\t". $f[4] ."\t".$desc."\t"."."."\t".$f[6]."\n";
	}
	elsif ($f_type eq 'interval'){
		my $desc = $f[0] .":".$f[1] ."-".$f[2].$f[3];
		print LIFT $f[0] ."\t". $f[1] ."\t". $f[2] ."\t".$desc."\t"."."."\t".$f[3]."\n";
	}
}
close(BED);
close(LIFT);