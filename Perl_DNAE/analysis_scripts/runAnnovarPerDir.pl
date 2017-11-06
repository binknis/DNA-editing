#Function: Runs annovar for a list of BED files and respective DBs found in two directories
#Input: 1. Directory of Annovar BED files (File suffix: .annovar.bed)
#		2. Dir of DBs for Annovar run
#		3. Dir for Annovar output
#		4. Gene database for annotation (gene (def)= Refseq, ensgene = Ensembl, knowngene = UCSC known genes)
#Note: will run for any bed file found in the BED-dir and present in the map list. 
use strict; 
(my $bedDir, my $dbDir, my $outDir, my $geneDB) = @ARGV; 
$geneDB = "gene" unless $geneDB; #default is Refseq; Use "ensgene" for Ensembl genes and "knowngene" for UCSC's known genes

##Get BED files 
opendir(BEDS, $bedDir) || print "Open Directory: $bedDir\n";
my @bedList = grep {/.annovar.bed/} readdir(BEDS);
closedir(BEDS);
#print "@bedList\n"; 

##Get assembly names
my $org2assem = "/home/alu/binknis/Analyses/Utils/Map_Clades_Orgs_Assembly_noHeader.txt"; #cols = Clade, Assembly, Organism.
my %o2a = (); 
open (MAP, $org2assem) or die "open $org2assem\n";
while (my $l = <MAP>){
	chomp $l; 
	(my $clade, my $org, my $assem) = split(/\t/, $l); 
	$o2a{$org} = $assem; 
}
close(MAP); 

##Run Annovar for each BED file
my $bedPre = "sites_A_"; 
my $bedSuff = "_LTR_bestParams.annovar.bed"; 
foreach my $org (sort keys %o2a){
	my $assembly = $o2a{$org}; 
	my $bedFile = $bedDir ."/". $bedPre .$org .$bedSuff; 
	if (-e $bedFile){ #BED file exists for this organism
		my $db = $dbDir ."/". $assembly ."db";
		system("/private/common/Software/ANNOVAR/annovar/annotate_variation.pl --geneanno --buildver " .$assembly ." -dbtype ". $geneDB ." ". $bedFile." ".$db);
	}
}
