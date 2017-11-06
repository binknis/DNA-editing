#Function: create sortGenomeN.sh scripts, by finding organisms in Raw_Data directory. 
#script1 must be run before script2 - for each run. 
use strict; 
my $rawDataDir = shift;
my @types; 
my $org; 
my $num=1; 
my $count = 0;

opendir(RAWDATA, $rawDataDir) || print "$rawDataDir didn't open\n";
@types = sort{lc($a) cmp lc($b)}(readdir(RAWDATA));
shift (@types) while ($types[0] =~ /^\./); #skip "." and ".." links
closedir(RAWDATA);

foreach my $type (@types){
	my %organisms=(); 
	my $typeDir = $rawDataDir ."/". $type; 
	next unless (-d $typeDir); #skip non-directory files
	
	opendir(TYPE, $typeDir) || print "$typeDir didn't open\n";
	my @orgList = readdir(TYPE); 
	closedir(TYPE);
	
	foreach my $file(@orgList){
		next if ($file =~ /^\./); #skip "." and ".." links
		next if (-d $file); #skip directory files
		if ($file =~ /^(\S+)\.(interval|fasta)$/){ 
			$organisms{$1} = 1; 
		}
	}

	open (UNMATCHED, ">>unmatchedOrganisms.txt") || die ("couldn't open unmatchedOrganisms.txt file\n"); 
	foreach $org (keys(%organisms)){
		my $intervalFile = $typeDir."/".$org.".interval"; 
		my $fastaFile = $typeDir ."/".$org.".fasta";
		if (-e $intervalFile and -e $fastaFile){
			my $script1 = "sortGenomes_".$$."_".$num.".sh"; 
			my $numPlusOne = $num+1; 
			my $script2 = "sortGenomes_".$$."_". $numPlusOne .".sh"; 
			open (SORT_SCRIPT1, ">>$script1") || print "file open $num failed\n";
			open (SORT_SCRIPT2, ">>$script2") || print "file open $num failed\n";
			#(my $seqDataFile, my $seqFile, my $organism, my @classes) = (shift, shift, shift, @ARGV);
			print SORT_SCRIPT1 "nohup perl sortGenome.pl $intervalFile $fastaFile $org LINE SINE LTR Other &\n"; 
			print SORT_SCRIPT2 "nohup perl sortGenome.pl $intervalFile $fastaFile $org &\n";
			close(SORT_SCRIPT1);
			close(SORT_SCRIPT2);
			chmod (0755, $script1, $script2);
			$count++; 
			$num += 2 if ($count % 10 == 0); 
		}
		else {
			print UNMATCHED $org . "\n";
		}
	}
	close(UNMATCHED); 
}
