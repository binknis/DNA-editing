#Function: calculate average divergence (mismatches, ins, del and total) for each subfamily in a file
#Input: 1. tableFile - created with for i in *.rmsk.txt.gz ; do  gunzip -c $i | sed -e 's/\?//g' -e 's/\//_/g' |  awk 'BEGIN {OFS="\t"} $12=="LTR" {print $12,$13,$11, $6 ":" $7 "-" $8 $10, $2,$3,$4,$5 }' > table_LTR_"$i"  ; done &
#					Which creates a table file for every rmsk.txt.gz file in 
#			Fields: [0]class, [1]fam, [2]subfam, [3]coords(chrN:\d+-\d+[+-] format), [4]swScore, [5]milliDiv, [6]milliDel, [7]milliIns
#Output: stdout with fields: [0]assembly, [1]class, [2]fam, [3]subfam, [4]ave_swScore, [5]ave_pcDiv, [6]ave_pcDel, [7]ave_pcIns, [8]ave_totalDiff, [9]num_sf_seqs
use strict; 
use Math::Round qw(nearest);

(my $tableFile, my $outDir) = @ARGV; 
(my $assembly) = $tableFile =~ /^table_LTR_(\S+).rmsk.txt.gz/; 
my %scores = (); #"0" - bitscore, "1" - milliDiv, "2" - milliDel, "3" - milliIns
my $colsPreScores = 4; 
my $scoreCols = 4; 

if ($tableFile =~ /\.gz$/){
	open (TABLE, "gunzip -c $tableFile |") or die "open $tableFile\n"; 
}
else{
	open (TABLE, $tableFile) or die "open $tableFile\n"; 
}
while (<TABLE>){
	chomp $_; 
	my @f = split(/\t/, $_); 
	my $taxa = join("\t", @f[0..2]); #class, fam, subfam
	#init if first time to see this taxa
	unless (exists $scores{$taxa}){
		my @arr = (0) x 4; 
		$scores{$taxa} = \@arr;
	}
	#add scores to sum
	for (0 .. ($scoreCols-1)){
		$scores{$taxa}[$_] += $f[$colsPreScores + $_]
	}
	$scores{$taxa}[$scoreCols]++; #count num of seqs for average
}
close(TABLE); 

foreach my $t (sort keys %scores){
	my $sumDivergence; 
	for my $pos(0 .. ($scoreCols-1)){ #divide by number of sequences
		$scores{$t}[$pos] = nearest(0.01, $scores{$t}[$pos] / $scores{$t}[$scoreCols] / 10); 
		$sumDivergence += $scores{$t}[$pos] if $pos; #sum all but SW score (pos=0)
	}
	print join ("\t", $assembly, $t, @{$scores{$t}}[0..3], $sumDivergence, $scores{$t}[4]), "\n";
}