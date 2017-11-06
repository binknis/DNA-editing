#function: For each organism in Raw_Data, calculate the amount of bps in various repeat classes
#		(currently calculates for DNA, LINE, LTR, SINE, but can be modified). 
#		Used to proceed and calculate the percent of each class in every genome for graphical presentation
#		(which is found in the "genomes" page in UCSC). 
#OUTPUT: print to file a tuple: Clade Organism Assembly	DNA	LINE LTR SINE  (tab delimited)
#	(Where DNA, LINE, LTR, SINE are number of nucleotides in each)
# Notes: 
#	1. need to make sure that all interval and fasta files in Raw_data are zipped
#	2. Run from Perl_DNAE file, one level under DNA_editing_results dir.

use strict;
use getAll; 
my $start; my $end; my $class; my $length;
my $raw_data = "../binknis_data/Raw_Data"; 
my @subdirs = qw(Deuterostome Insect Mammalian Nematode Other Vertebrate); 
my @classList = qw(DNA LINE LTR SINE);
my $classRegex = "(".join('|',@classList).")";
open (ERRS, ">../DNA_editing_results/calcAmountOfRT_errors.txt") || die "open error file failed\n"; 
open (OUT, ">../DNA_editing_results/amount_RT_bp_DNA_LINE_LTR_SINE.txt") || die "open out file failed\n"; 
foreach my $clade (@subdirs){
	my $cladeDir = $raw_data ."/" .$clade; 
	opendir(CLADE, $cladeDir) || print "$cladeDir didn't open\n";
	my @fileList = sort{lc($a) cmp lc($b)} grep {/\.interval.gz$/} (readdir(CLADE));
	shift(@fileList) while ($fileList[0] =~ /^\./); #erase "." and ".." links
	closedir(CLADE);	
	#calculate the amount of bp for each class in each organism
	foreach my $file (@fileList){
		#get and calculate the data
		my %totalLen = (); 
		(my $org, my $assembly) = getOrgAndAssembly($cladeDir, $file); 
		if ($org eq "" && $assembly eq ""){
			print ERRS "$cladeDir/$file doesn't have a matching fasta!\n"; 
			next; 
		}
		
		my $intervalFile = $cladeDir."/".$file; 
		open (INTERVAL, "gunzip -c $intervalFile |") || die "$intervalFile didn't open\n"; 
		while (my $line = <INTERVAL>){
			next if $line =~ /^#/;
			# print "line: $line\n"; 
			chomp $line;
			# print "he\n"; 
			my @data = split (/\s+/, $line); 
			($start, $end, $class) = ($data[1], $data[2], $data[5]); 
			next unless $class =~ /$classRegex/; 
			$length = $end - $start; 
			$totalLen{$class}=0 unless exists $totalLen{$class}; 
			$totalLen{$class} += $length; 
		}
		close(INTERVAL);
		
		# print the output
		print OUT $clade."\t".$org."\t".$assembly; 
		foreach my $class (@classList){
			print OUT "\t0" unless exists $totalLen{$class}; 
			print OUT "\t".$totalLen{$class}; 
		}
		print OUT "\n"; 
	}
}
close (ERRS);
close (OUT);  
#get the organism name (from the file name) and the assembly name (from the first defline in the fasta file)
sub getOrgAndAssembly(){
	(my $cladeDir, my $interval) = @_;
	my $assembly; 
	$interval =~ s/\.interval.gz$//g; 
	my $org = $interval; 
	my $fasta = $cladeDir ."/". $org . ".fa.gz"; 
	open (FASTA, "gunzip -c $fasta |") || return ("",""); 
	my $line; 
	while ($line = <FASTA>){
		if ($line =~ /^>/){
			($assembly) = $line =~ />([^_]+)/; 
			# print $org ."\t". $assembly."\n"; #optional if you want to print the name mapping
			last; 
		}
	}
	close(FASTA); 
	return ($org, $assembly); 
}
