#Function: runs intersect bed between beds for specified assemblies with its respective BED files (refseq, ensgene - mRNA, CDS, intron, exon, codingExon etc.)
(my $bedDir, my $alldbdir, my $assemblies) = @ARGV; 

$assemblies =~ s/,$//; #enable trailing comma in input
my @assembly_arr = split (',', $assemblies); 
my $assembly_regex = join('|', $assembly_arr); 
$assembly_regex = "^(" . $assembly_regex . ")$";
print $assembly_regex ."\n"; 

foreach my $as (){
	my $a_file = $bedDir; 
	foreach 
	$b_file = ; 

	


}
