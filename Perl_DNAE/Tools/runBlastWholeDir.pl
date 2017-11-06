use strict; 

(my $queryDir, my $dbFile, my $outFile, my $blastType, my $cores, my $pval, my $S) = @ARGV; 

opendir (DNA, $queryDir) || die "couldn't open $queryDir\n"; 
my @queryFiles = sort{lc($a) cmp lc($b)}(readdir(DNA)); 
shift(@queryFiles) while ($queryFiles[0] =~ /^\./);
close (DNA); 

foreach my $queryFile(@queryFiles){
	#(my $DNAFile, my $dbFile, my $outFile, my $cores) = @ARGV; 
	my $runBlast = "perl runBlast.pl $queryDir/$queryFile $dbFile $outFile $blastType $cores $pval $S";
	$runBlast =~ s/\(/\\\(/g; $runBlast =~ s/\)/\\\)/g; #add backslashes before "(" and ")" symbols.
	system ($runBlast);
}