#output is appended to output file
#NOTE: input pval must have '1e-' prefix (unless you wanted to enter an integer for pval)
#		2. if outFile eq stdout will print to stdout. 
#		3. -K parameter was disabled
use strict; 

(my $queryFile, my $dbFile, my $outFile, my $blastType, my $cores, my $pval, my $S, my $m) = @ARGV; #S is strand, m is outfmt
my $p;  
my $formatdb; 
my $blast; 
# $pval = "1e-" . $pval unless $pval =~ /1e-/;

if ($blastType eq 'blastn'){
	$pval = "1e-50" unless ($pval); #DEFAULT PVAL
	$p = 'F'; 
}
elsif ($blastType eq 'blastx'){
	$pval = "1e-20" unless ($pval); #DEFAULT PVAL
	$p = 'T'; 
}
elsif ($blastType eq 'blastp'){
	$pval = "1e-20" unless ($pval); #DEFAULT PVAL
	$p = 'T'; 
}
else {die "undefined blastType: $blastType\n";}

$cores = 1 if ($cores == 0); 
$S = 3 if ($S == 0);

# system("echo $pval"); 
#mkdir "blastOutDir"; 
#$queryFile =~ /([^\/]+)$/; 
#$outFile  = "blastOutDir/$outFile" . "_" . $1; #***

$formatdb = "formatdb -i $dbFile -p $p -o T\n" unless -e $dbFile .".nhr";
$blast = "blastall -p $blastType -d $dbFile -i $queryFile -e $pval -F F -a $cores -S $S -v 0 "; 
# $blast .= "-K $K " if ($K); 
$blast .= "-m $m " if ($m); 
$blast .= ">> $outFile" unless ($outFile eq 'stdout'); 

$formatdb =~ s/\(/\\\(/g; $formatdb =~ s/\)/\\\)/g; #add backslashes before "(" and ")" symbols. 
$blast =~ s/\(/\\\(/g; $blast =~ s/\)/\\\)/g; #add backslashes before "(" and ")" symbols. 

system ($formatdb); 
system ($blast);
