#Function: get the BED files of retroelements based on a list of organisms, class, subfamily and (optionally) subfamilies
#Note: the input BED file is in the following format (contains taxa in NAME field)
#   	scaffold_0      364     478     LTR|ERVL-MaLR|MLT1B     .       -
##Example input files and dir: 

#$famListFile --- /home/alu/binknis/Rscripts/DNAE/Choose_best_params/lowConfidence_byClusters_AmoreDiv_mostGMapGinPairs/Tracks/orgClassFam_lowConfidence_byClusters_AmoreDiv_mostGMapGinPairs.txt
#$famListFile --- /home/alu/binknis/Rscripts/DNAE/Choose_best_params/lowConfidence_byClusters_AmoreDiv_mostGMapGinPairs/Tracks/orgClassFamSubfam_lowConfidence_byClusters_AmoreDiv_mostGMapGinPairs.txt
#$bedDir --- /private2/Dropbox/users/binknis_data/Raw_Data/BEDfiles/LTR

use strict; 
(my $famListFile, my $bedDir, my $outDir) = @ARGV; 

my $subfamsFound=0; 
#Get list of fams or subfams
my %taxas = (); 
open (my $fams_fh, $famListFile) or die "open $famListFile\n"; 
while(my $l = <$fams_fh>){
	chomp $l; 
	(my $org, my $class, my $fam, my $subfam) = split("\t", $l ); 	my $taxa = ($subfam ? join('|', $class, $fam, $subfam) : join('|', $class, $fam)); 
	$taxas{$org}{$taxa} = 0; 
	$subfamsFound++ if ($subfam); 
}
close($fams_fh); 

#Get rows from BED file based on taxa (if exists in specified list of fams/subfams)
foreach my $org (sort keys %taxas){
	my $orgBed_file = $bedDir ."/". $org .".bed"; 
	my $orgOut_file = $outDir ."/". $org .".bed"; 
	open (OUT, ">".$orgOut_file) or die "open $orgOut_file\n"; 
	open (BED, $orgBed_file) or die "open $orgBed_file\n"; 
	while (my $l = <BED>){
		chomp $l; 
		my @f = split("\t", $l); 
		my $taxa = $f[3]; 
		$taxa =~ s/\|[^|]+$// unless($subfamsFound); 
		if (exists $taxas{$org}{$taxa}){
			print OUT $l ."\n"; 
		}
	}
	close(OUT); 
	close(BED); 
}
