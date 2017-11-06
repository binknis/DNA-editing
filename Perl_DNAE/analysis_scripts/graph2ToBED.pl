#Function: create a bed file from the graph2 file (was written to create a BED file with the taxa of the element)

use strict; 
(my $graph2_file, my $mm, my $SorT) = @ARGV; 
$mm = "GA" unless $mm;
my $bed_file = $graph2_file; 
my $replacement = "seq_" . $editNuc; 
$bed_file =~ s/graph2/$replacement/;
open (my $bed_fh, ">" . $bed_file) or die "open $bed_file\n"; 
while (my $l = <$g_fh>){
	chomp $l; 
close($bed_fh); 