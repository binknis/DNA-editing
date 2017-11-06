#Function: creates a file with nucleotide frequencies for a specific family
#Output: file - "../Data/$org/$class/db/NucByFams.txt"; one line per family. 
use strict; 
use getAll; 

(my $org, my $class) = @ARGV; 
#open and truncate existing file
my $fam_nuc_freq_file = "../Data/$org/$class/db/NucByFams.txt";
open (OUT, ">$fam_nuc_freq_file") || die "couldn't open $fam_nuc_freq_file\n";

my $fams = getAll::families($org, $class); 
for my $fam(@$fams){
	#get num bps by name
	my $elemAndBPDir = "../DNA_editing_results/db_stats/ElementsAndBasepairs_allOrgs"; 
	my $name_elemBP_file = $elemAndBPDir ."/allOrgs_bySubfamily.txt";
	my $name_elemBP = getAll::lines($name_elemBP_file); 
	my $regex = $org ."\t". $class ."\t". $fam;
	my @name_counts = grep (/$regex/ , @$name_elemBP); 

	#calculate a weighted mean by:
	#1. per name, multiply freq by num seqs and summing 
	# my @nucs = qw/a c g t/; 
	my @nucs_freq = (0) x 4;
	for (@name_counts){
		(my $subfam, my $bps) = $_ =~ /(\S+)\t\d+\t(\d+)\s*$/;
		my $nuc_subfam = getAll::nucFreqSubfam($org, $class, $fam, $subfam); #get subfam's Nuc freq
		for my $i(0 .. 3){
			$nucs_freq[$i] += $bps * $nuc_subfam->[$i];
			# print $subfam ." ". $nuc_subfam->[$i] ."\n"; 
		}
	}

	#2. divide by num seqs in fam
	my $fam_elemBP_file = $elemAndBPDir ."/allOrgs_byFamily.txt";
	my $fam_elemBP = getAll::lines($fam_elemBP_file); 
	my $regex = $org ."\t". $class ."\t". $fam; 
	(my $my_fam_seq_count) = grep (/$regex/ , @$fam_elemBP);
	(my $bps) = $my_fam_seq_count =~ /(\d+)\s*$/;
	for my $i(0 .. 3){
		if ($bps > 0){
			$nucs_freq[$i] /= $bps; 
		}
		else{
			print "no bps found in $fam_elemBP_file file for $org $class $fam\n"; 
			next; 
		}
	}
	
	#print output
	print OUT $org ."\t". $class ."\t". $fam ."\t"; 
	# print $org ."\t". $class ."\t". $fam ."\t"; 
	for my $i(0 .. 3){
		print OUT $nucs_freq[$i] . ($i == 3 ? "\n" : "\t"); 
		# print $nucs_freq[$i] . ($i == 3 ? "\n" : "\t"); #***
	}
}
close(OUT);