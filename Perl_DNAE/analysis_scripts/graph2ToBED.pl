#Function: create a bed file from the graph2 file (was written to create a BED file with the taxa of the element)

use strict; 
(my $graph2_file, my $mm, my $SorT) = @ARGV; 
$mm = "GA" unless $mm;my @mmAr = split('', $mm);  my $editNuc = $SorT ? $mmAr[0] : $mmAr[1]; #default target, if 1 is specified then for source
my $bed_file = $graph2_file; 
my $replacement = "seq_" . $editNuc; 
$bed_file =~ s/graph2/$replacement/;$bed_file =~ s/\.txt$/.bed/;open (my $g_fh, $graph2_file) or die "open $graph2_file\n"; 
open (my $bed_fh, ">" . $bed_file) or die "open $bed_file\n"; 
while (my $l = <$g_fh>){
	chomp $l; 	(my $org, my $class, my $fam, my $subfam, my $seqS, my $seqT) = my @f = split("\t", $l); 	my $wantedCoords = $SorT ? $seqS : $seqT;	$wantedCoords =~ /^([^:]+):(\d+)-(\d+)([+-])$/;	my $taxa = join('|', $class, $fam, $subfam);	print $bed_fh join("\t", $1, $2, $3, $taxa, ".", $4) . "\n";}
close($bed_fh); close($g_fh); 