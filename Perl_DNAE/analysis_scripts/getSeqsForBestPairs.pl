#Function: gets sequences of editing-pairs from the respective sequence files
use strict; 

(my $clustFile, my $outDir) = @ARGV; 

my $home =  $ENV{"HOME"}; 

($outDir) = $clustFile =~ /(.*)\// unless $outDir; 
$outDir .= "/editingPairSequences"; 
mkdir $outDir; 

#Get Assembly to Organism map.
#fields: clade org assembly
my $assemblyToOrgFile = "/home/alu/binknis/Analyses/Utils/Map_Clades_Orgs_Assembly.txt";
my %a2o = ();
open (ATO, $assemblyToOrgFile) or die "$assemblyToOrgFile\n";
while (my $l = <ATO>){
        chomp $l;
        my @f = split(/\t/, $l);
        $a2o{$f[2]} = $f[1];
}
close(ATO);

#Read data from clusters and fetch sequences for each coordsS and coordsT pair
open (CLUSTS, $clustFile) or die "open $clustFile\n"; 
while(my $l = <CLUSTS>){
	next if $l =~ /^mismatch/; 
	#get pair info
	chomp $l; 
	my @f = split(/\t/, $l); 
	(my $mm, my $assembly, my $class, my $fam, my $sf, my $coordsS, my $coordsT) = ($f[0], $f[1], $f[2], $f[3], $f[4], $f[5], $f[6]); 
	
	#open db seq file and fetch wanted sequences
	my @regexes = ($coordsS, $coordsT); 
	my $outFile = $outDir ."/". join("_", $mm, $a2o{$assembly}, $class, $fam, $sf, $coordsS, $coordsT); 
	$outFile =~ s/:/_/g; 
	my $dbSeqFile = join ("/", $home, "Data", $a2o{$assembly}, $class, "db", "files_".$fam, "Seq_".$sf); 
	open (OUT, ">>" . $outFile) || die "open $outFile\n";
	open (FASTA, $dbSeqFile) || die "couldn't open $dbSeqFile\n";
	my $seq;
	while (my $line = <FASTA>){
		if ($line =~ /^>/){
			$seq = $line;
			while ($line = <FASTA>){
				if ($line =~ /^>/){
					seek(FASTA, -length($line), 1);
					last;
				}
				$seq .= $line;
			}
		}
		(my $seqID) = ($seq =~ /^>(.+)?\n/); #get seqs ID from head of defline
		chomp $seq;
		# print $seqID ."\n";
		#print the sequence only if it matches one of the regexes
		my $found = 0;
		foreach my $regex(@regexes){
			$found = 1 if  $seqID =~ /$regex/;
		}
		(my $defline, my $seqLines) = split("\n", $seq, 2); 
		print OUT $defline . "\n" . lc($seqLines) . "\n" if $found;
	}
	close (FASTA);	
	close (OUT);
}
close(CLUSTS); 
