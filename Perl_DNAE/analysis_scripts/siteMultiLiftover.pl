#Function: liftover set of intervals to list of other genomes
#Example for assembly list: CalJac3,GorGor3,MicMur1,NomLeu1,OtoGar1,PanTro3,PonAbe2,RheMac2,RheMac3,TarSyr1
#Input: 1. $bedFile - bed or gff (if isn't in good format for liftover, where the name is 
#		2. assembly: to lift from (typically, list of edited sites)
#		3. newAssemblies: comma-separated list of editing sites to lift to
use strict; 

## consts
my $home = $ENV{HOME}; 
my $genomeDir = "/home/alu/binknis/binknis_data/DBs/UCSC/genomes"; 
my $chainDir = "/home/alu/binknis/binknis_data/DBs/UCSC/liftOver/over.chains"; 
#my $outFile; see below
## args
(my $bedFile, my $assembly, my $newAssemblies, my $outFile) = @ARGV; 
my @newAs = split (',', $newAssemblies); 

## validate that input bed file is good for liftover input (has coords in name so that you see the original coords in lifted output file)
unless ($bedFile =~ /forLift\.bed/){
	system("perl516 ". $home ."/Perl_DNAE/Tools/convertToLiftoverBed.pl $bedFile"); 
	$bedFile =~ s/(bed|gff)$/forLift.bed/;
}
unless ($outFile){
	my $outFile = $bedFile;
	$outFile =~ s/forLift\.bed/nucLiftRes.tab/; 
}
## liftover and get nucleotides in respective positions
#get all coords from original file
my $coord; 
my $coords;
my @coordsList = (); #save order of editing sites, for printing output in same order
open (COORDS, $bedFile) || die "open $bedFile\n"; 
while (my $l = <COORDS>){
	chomp $l; 
	my @fields = split (/\t/,$l); 
	push (@coordsList, $fields[3]); 
}
close(COORDS); 
#init hash for all assemblies - zero for each assembly (will be replaced for successfully lifted sites)
my %nucFound = ();
foreach my $na (@newAs){	
	foreach my $c (@coordsList){
		$nucFound{$na}{$c}=0; 
	}
}

#get nucs
foreach my $na (@newAs){	
	## liftover
	#example command: "liftOver input.bed hg18ToHg19.over.chain.gz output.bed unlifted.bed"
	(my $firstChar, my $rest) = $na =~ /(\S)(\S+)/;
	my $fromTo = $assembly ."To". (uc $firstChar) .$rest; 
	my $chainFile = $chainDir ."/". $fromTo .".over.chain.gz" ; 
	my $lifted = $bedFile .".".$fromTo.".lifted"; 
	my $unlifted = $bedFile .".".$fromTo.".unlifted";
	my $lift = "/home/alu/binknis/binknis_data/DBs/UCSC/liftOver/liftOver $bedFile $chainFile $lifted $unlifted";
	system($lift); 
	
	##map old assembly coords to new assembly coords (= lift results to lift input)
	my %coordMap = (); 
	open (LIFTED, $lifted) || die "open $lifted\n"; 
	while (my $l = <LIFTED>){
		chomp $l; 
		my @f = split (/\t/,$l); 
		my $new = $f[0] .":". $f[1] ."-". $f[2] . $f[5]; 
		my $old = $f[3]; 
		$coordMap{$new} = $old;
		# print "old " . $old ." new ". $new ."\n"; 
	}
	close(LIFTED);
	
	## get nucs
	my $genomeFile = $genomeDir ."/". $na .".fa"; 
	my $liftedTab = $lifted .".tab"; 
	my $getFastas = "bedtools getfasta -s -tab -fi $genomeFile -bed $lifted -fo $liftedTab"; 
	system ($getFastas); 
	
	## insert lifted nucs to respective hash
	#Note: addressed in output parsing: a. parenthases b. nucs to lc 
	open (LIFTNUC, $liftedTab) || die "open $liftedTab\n"; 
	while (my $l = <LIFTNUC>){
		chomp $l; 
		$l =~ s/[()]//g; #erase parentheses added to coords (around strand sign) in bedtools getfasta output 
		(my $coord, my $nuc) = $l =~ /^(\S+)\t(\w+)/;
		# print $coord ."\t". $nuc ."\n"; 
		#print  $coord ."\t". $coordMap{$coord} ."\n"; #***
		$nucFound{$na}{$coordMap{$coord}} = lc $nuc; 
	}
	close(LIFTNUC);
	
	#delete intermediate files
	unlink($lifted, $unlifted, $liftedTab); 
}

#print results in tabular format to file
open (OUT, ">$outFile") or die "open $outFile\n"; 
print OUT "#Coords\t" , join('\t', @newAs), "\n"; 
foreach my $c (@coordsList){
	print OUT $c; 
	foreach my $na (@newAs){
		print OUT "\t" . $nucFound{$na}{$c}; 
	}
	print OUT "\n";
}
close(OUT); 