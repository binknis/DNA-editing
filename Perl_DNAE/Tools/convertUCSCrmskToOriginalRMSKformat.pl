#Function: converts UCSC's table browser's RMSK table to fit the original RMSK output format

### The formats and needed chanegs
### UCSC rmsk columns: bin, swScore, milliDiv, milliDel, milliIns, genoName, genoStart, genoEnd, genoLeft, strand, repName, repClass, repFamily, repStart, repEnd, repLeft, id
#Example:
#607     12955   105     9       10      chr1    3000000 3002128 -192469843      -       L1_Mus3 LINE    L1      -3055   3592    1466    1
### Original rmsk columns ([*] = format change needed for column): swScore, pcDiv, pcDel, pcIns, genoName, genoStart, genoEnd, genoLeft[*], strand, repName, repClass/repFamily, repStart[*], repEnd, repLeft[*], id
#Headers:
# SW  perc perc perc  query    position in query     matching repeat      position in  repeat
# score div. del. ins.  sequence begin  end  (left)    repeat  class/family   begin end (left) ID
#Example:
#1320 15.6  6.2  0.0  HSU08988  6563 6781 (22462) C  MER7A   DNA/MER2_type    (0)  337  104  20

### List of format changes needed: 
#1. delete bin column
#2. change format of genLeft (instead of negative value make abs and add parentheses) 
#3. Same as (2) for repStart or repLeft, when negative

use strict; 

(my $rmskFile, my $outFile) = @ARGV; 

if($rmskFile =~ /\.gz$/){
	open(RMSK, "gunzip -c $rmskFile |") or die "can't pipe $rmskFile\n"; 
}
else{
	open(RMSK, $rmskFile) or die "open $rmskFile\n"; 
}

open(my $rmskOut, ">" . $outFile) or die "open $outFile\n"; 
while(my $l = <RMSK>){
	chomp $l; 
	next if $l =~ /^\s*#/; #skip comment lines
	next if $l =~ /^\s*$/; #skip empty lines
	next if $l =~ /^\s*(SW|score)/; #skip header lines
	#UCSC rmsk columns: bin, swScore, milliDiv, milliDel, milliIns, genoName, genoStart, genoEnd, genoLeft, strand, repName, repClass, repFamily, repStart, repEnd, repLeft, id
	my @fs = split("\t", $l); 
	die "bad number of columns: $#fs . Should be 17 for UCSC tables\n" unless ($#fs == 16); #has bin column (this should happen...)
	(my $bin, my $swScore, my $milliDiv, my $milliDel, 
		my $milliIns, my $genoName, my $genoStart, my $genoEnd, my $genoLeft, my $strand, 
		my $repName, my $repClass, my $repFamily, my $repStart, my $repEnd, my $repLeft, my $id) = @fs;
	
	#convert milli to percent
	my $pcDiv = $milliDiv / 10; 
	my $pcDel = $milliDel / 10; 
	my $pcIns = $milliIns / 10; 
	
	#change format of toEnd columns
	if($genoLeft < 0){
		$genoLeft = "(" . abs($genoLeft) .")"; 
	}
	if($repStart < 0){
		$repStart = "(" . abs($repStart) .")"; 
	}
	if($repLeft < 0){
		$repLeft = "(" . abs($repLeft) .")"; 
	}
	
	#change strand
	if($strand eq "-"){
		$strand = "C"; 
	}
	
	my $repClassFam = $repClass ."/". $repFamily; 
	
	my $outLine = join("\t", $swScore, $pcDiv, $pcDel, $pcIns, $genoName, $genoStart, $genoEnd, $genoLeft, $strand, $repName, $repClassFam, $repStart, $repEnd, $repLeft, $id);
	print $rmskOut $outLine ."\n";
}
close(RMSK); 
close($rmskOut);
