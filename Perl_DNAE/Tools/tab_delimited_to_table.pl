#converts a cluster_stats file to tabular format (Threshold X Pvalue: for any other field)
#Th is rows, Pval is columns

use strict; 
my $line; 
my %hash; 
(my $file, my $field_num) = @ARGV; 
open (CSTATS, "<$file") || die "open $file\n"; 

while ($line = <CSTATS>){
	chomp $line; 
	next if $line !~ /\S+/; 
	my @fields = split (/\t/, $line); 
	$hash{$fields[1]}{$fields[2]} = $fields[$field_num]; 
}
close(CSTATS); 
 
foreach my $pval(sort{$b <=> $a} keys %hash){
	print "\t$pval"; 
}

my @pvals = sort{$b <=> $a} keys %hash; 
my @thresholds = sort{$a <=> $b} keys %{$hash{$pvals[0]}}; 

foreach my $th (@thresholds){
	print "\n".$th; 
	foreach my $pval (@pvals){
		print "\t". $hash{$pval}{$th}; 
	}
}

print "\n"; 