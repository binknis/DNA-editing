use Bio::SearchIO;
#FUNCTION: Written to get sequences from a set of sequences blasted against a reference genome. 
#			Was designed to get L2 sequences from the Anolis genome, using 3 reference sequences. 
#INPUT: 1. Blast of sequences against the genome. 
# 		2-3. Interval and sequence file names - full path (for output)
#OUTPUT: creates two files in the RAW_DATA format: interval, fasta. 

use Bio::Root::Root;
use List::Util qw[min max];

(my $blastFile, my $intervalFile, my $seqFile) = @ARGV; 
open (INVL,">$intervalFile") || die "didn't open $intervalFile\n"; 
open (SEQ,">$seqFile") || die "didn't open $seqFile\n"; 

my $genome = "anoCar2"; 

my $in = new Bio::SearchIO(-format => 'blast', -file   => $blastFile);
my $counter = 0; 
my $seq; 
while( my $result = $in->next_result )
{
	while( my $hit = $result->next_hit )
	{
		while( my $hsp = $hit->next_hsp )
		{	
			# print ">" . $hit->name . ":" . $hsp->start('hit') . "-" . $hsp->end('hit'); #write defline (position as "chr:start-end")
			print INVL $hit->name ."\t". $hsp->start('hit') ."\t". $hsp->end('hit') ."\t". (($hsp->strand('hit') == 1) ? "+" : "-")."\t". $result->query_name ."\t". "LINE" ."\t". "L2" . "\n";  
			print SEQ ">" . $genome ."_". $hit->name . "_" . $hsp->start('hit') . "_" . $hsp->end('hit') ."_". (($hsp->strand('hit') == 1) ? "+" : "-") . "\n"; #write defline (position as "chr:start-end[+-]")
			$seq = $hsp->hit_string; 
			$seq =~ s/-//g; #remove gaps
			my $pos = 0; 
			my $length = length($seq); 
			while ($pos < $length){
				print SEQ substr ($seq, $pos, min(80,$length-$pos)). "\n";
				$pos += 80; 
			}
			# print $hsp->hit_string . "\n"; 
			# exit if $counter++ == 10; 
		}
	}
}

close(INVL);
close(SEQ); 