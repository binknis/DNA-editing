#Function: create frequency table from sequences of same length.
#Input: 'edited_pos': the edited position in the sequence (1-base!); 
#		'margin': how many nucs upstream and downstream to fetch.
#		normFlag: 0 (default) - return fraction, 1 - return frequency 
#Note: Logo file must have ".txt" suffix; Output will have same name with "_freq.txt" suffix
#MAKE SURE THIS WORKS - I ADDED LINES 44-46 (which avoid division by zero). 
use strict; 
(my $logo_file, my $edited_pos, my $margin, my $normFlag) = @ARGV; 
my %hist = (); 
my @nucs = ('a','c','g','t');
$margin = $edited_pos - 1 unless $margin;
my $num_poses = 2 * $margin + 1; 
$edited_pos--; #-1 for 0-base. 
my $start = $edited_pos - $margin; 
my $end = $edited_pos + $margin;
my $seq_count=0; 

for (my $i=$start; $i<=$end; $i++){
	foreach my $nuc(@nucs){
		$hist{$i}{$nuc}=0; 
	}
}

open (LOGO, "<".$logo_file) || die "couldn't open logo file\n"; 
while (my $line = <LOGO>){
	chomp $line; 
	next if $line eq ''; #skip empty lines; 
	$seq_count++; 
	$line = lc $line;
	my @str = split('',$line);
	for (my $i=$start; $i<=$end; $i++){
		$hist{$i}{$str[$i]}++;
	}
}
close(LOGO); 

my $logo_freq_file = $logo_file; 
$logo_freq_file =~ s/\.txt$//;
if ($normFlag){
	$logo_freq_file .= '_freq2.txt'; 
}
else{
	$logo_freq_file .= '_freq.txt'; 
}

open (LOGO_FREQ, ">$logo_freq_file") or die "couldn't open $logo_freq_file\n"; 
print  LOGO_FREQ "position\tA\tC\tG\tT\n";
foreach my $pos (sort {$a <=> $b} keys %hist){
	print LOGO_FREQ ($pos - $edited_pos) . "\t";
	foreach my $nuc (@nucs){
		if ($seq_count == 0){
			$hist{$pos}{$nuc}=0;
		}
		else{
			$hist{$pos}{$nuc} /= $seq_count unless $normFlag; 
		}
		print  LOGO_FREQ $hist{$pos}{$nuc} . ($nuc eq 't' ? "\n" : "\t");
	}
}
close(LOGO_FREQ);