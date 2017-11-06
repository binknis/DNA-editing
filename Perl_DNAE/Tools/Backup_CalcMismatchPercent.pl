use strict; 
use IO::Compress::Gzip qw(gzip);
use IO::Uncompress::Gunzip qw(gunzip);


my $organism = @ARGV; 
my %mmPercent = (); #will be a hash of array 
my %mmPercentClass = (); 

my $classDir = "../Data/" . $organism; 
opendir(CLASSES, $classDir) || print "classdir didn't open\n";
my @classList = sort{lc($a) cmp lc($b)}(readdir(CLASSES));
shift(@classList) while ($classList[0] =~ /^\./); #erase "." and ".." links
closedir(CLASSES);

foreach my $class (@classList) {
	my $familyDir = $classDir ."/$class/results/blasts"; 
	opendir(FAMILY, $familyDir) || print "Family directory for $familyDir didn't open\n";
	my @familyList = sort{lc($a) cmp lc($b)}(readdir(FAMILY));
	shift(@familyList) while ($familyList[0] =~ /^\./); 
	closedir(FAMILY);
	foreach my $family (@familyList){
		my $nameDir = $familyDir . "/$family";
		opendir(NAME, $nameDir);
		my @nameList =  sort{lc($a) cmp lc($b)}(readdir(NAME));
		shift(@nameList) while ($nameList[0] =~ /^\./); 
		closedir(NAME);
		foreach my $name (@nameList) {
			getNamesMMPercent(); 
		}
	}
	printFamMMPercent();
}

sub getNamesMMPercent{
	#unzip and parse blast file
	my $blast_archive = $nameDir."/".$name . '.gz';
	my $blast_fileName = "blastForMMCalc".$$.".txt";
	unless (gunzip $blast_archive => $blast_fileName){
		print "gunzip failed for $blast_archive\n";
		return; 
	}
	#get mm data
	open(BLAST, $blast_fileName) || print "$blast_fileName didn't open!\n"; 
	while my $line (<BLAST>){
		if ($line =~ /^\s+Identities = \d+\/\d+ \((\d+)%\)/){
			my $mmPercent = 100 - $1;
			%mmPercent{$family}{$mmPercent}=0 unless exists %mmPercent{$family}{$mmPercent};
			%mmPercent{$family}{$mmPercent}++;
			%mmPercentClass{$mmPercent}=0 unless exists %mmPercentClass{$mmPercent};
			%mmPercentClass{$mmPercent}++;
		}
	}
	close(BLAST); 
	unlink($blast_fileName);
}

#print Classe's mm data to its file
sub printFamMMPercent{
	my $mmPercentFile = ">>../Data/" . $organism . "/" . $class . "/results/mmPercent.txt";
	open (OUT, $mmPercentFile) || print "$mmPercentFile didn't open\n"
	print OUT "Total for $class: \n";
	foreach my $percent(sort{$a <=> $b}(keys(%mmPercentClass))){
		print OUT $percent . "=" . $mmPercentClass{$percent} . " ";
	} 
	print OUT "\nHistogram by families: \n"; 	
	print OUT "$family:\n"; 
	
	foreach my $family (@familyList){
		foreach my $percent(sort{$a <=> $b}(keys(%mmPercent{$family}))){
			print OUT $percent . "=" . $mmPercentClass{$percent} . " ";
		}	
	}
	close(OUT); 
}