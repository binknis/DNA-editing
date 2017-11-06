#Function: converts a clusters file to tabular format (tuple per cluster)
#Input: cluster file (full-path), tabular file (full path)
#Output: creates the specified tabular file.
use strict;
(my $cluster_file, my $tabular_file) = @ARGV;
open (CLUSTS ,"<" . $cluster_file) || die ("couldn't open $cluster_file"); 
open (TAB ,">" . $tabular_file) || die ("couldn't open $tabular_file"); 

my $next_line;

my @arrayA = ();
my @arrayG = ();

#--table------
my $assembly; 
my $class;
my $family;
my $name;

my $chrG;
my $strandG;
my $startG;
my $endG;

my $chrA;
my $startA;
my $endA;
my $strandA;

my $mmSerials; 
my $locG="";
my $locA="";

my $align_length;
my $direct_mms;
my $all_mms; 
my $prob; 
my $num_clusts; 

#--------------
my $tuple;
while(my $line  = <CLUSTS>)
{
	if($line =~ /^Found/)
	{
		### READ CLUSTER ### 
		my $g_name_l  = <CLUSTS>;
		chomp $g_name_l; 
		my $a_name_l = <CLUSTS>;
		chomp $a_name_l; 
		my $prob_clustNum_l = <CLUSTS>;
		chomp $prob_clustNum_l; 
		my $mm_serial_num_l = <CLUSTS>;
		chomp $mm_serial_num_l; 
		my $locG_l = <CLUSTS>;
		chomp $locG_l; 
		my $locA_l = <CLUSTS>;
		chomp $locA_l; 
		my $len_mmNum_l = <CLUSTS>;
		chomp $len_mmNum_l; 
	
		### G name line ###
		#G name = 1082=hg19=chr1:144836998-144837209-=LINE=L1=HAL1
		($assembly, $chrG, $startG, $endG, $strandG, $class, $family, $name) = $g_name_l =~ /G name = \d+=([^=]+)=([^:]+):(\d+)-(\d+)([+-])=([^=]+)=([^=]+)=(\S+)/; 
	
		### A name line ###
		# A name = 118=hg19=chr1:16882278-16882479+=LINE=L1=HAL1
		($assembly, $chrA, $startA, $endA, $strandA, $class, $family, $name) = $a_name_l =~ /A name = \d+=([^=]+)=([^:]+):(\d+)-(\d+)([+-])=([^=]+)=([^=]+)=(\S+)/;
		
		### total_prob line ###
		# Total probability = 0.000172152888096637, 1 cluster
		($prob, $num_clusts) = $prob_clustNum_l =~ /Total probability = (\S+), (\d+) cluster/; 		
		
		### Edited MM serial line ###
		# Edited MM serial no.    :     1     2     3     4     5
		($mmSerials) = $mm_serial_num_l =~ /Edited MM serial no.    :\s+(.+\d+)/; 
		$mmSerials =~ s/\s+/,/g;
		
		### G locations line ###
		# Locations G seq:    72    77    89    99   141
		($locG) = $locG_l =~ /Locations G seq:\s+(.+\d+)/; 
		$locG =~ s/\s+/,/g;
		
		### A locations line ### 
		# Locations A seq:    71    76    88    98   140
		($locA) = $locA_l =~ /Locations A seq:\s+(.+\d+)/; 
		$locA =~ s/\s+/,/g;
		
		### Length, Mismatches ### 
		# Total length: 201, Direct mismatches: 5, All mismatches: 16
		($align_length, $direct_mms, $all_mms) = $len_mmNum_l =~ /Total length: (\d+), Direct mismatches: (\d+), All mismatches: (\d+)/;
	
		my $tuple = $assembly."\t"
					.$class."\t"
				   .$family."\t"
				   .$name."\t"
				   .$chrG."\t"
				   .$startG."\t"
				   .$endG."\t"
				   .$strandG."\t"
				   .$locG."\t"
				   .$chrA."\t"
				   .$startA."\t"
				   .$endA."\t"
				   .$strandA."\t"
				   .$locA."\t"
				   .$align_length."\t"
				   .$direct_mms."\t"
				   .$all_mms."\t"
				   .$mmSerials."\t"
				   .$prob."\t"
				   .$num_clusts; 
		
		# print $data "$tuple\n";
		print TAB "$tuple\n";
	}
}
 
close (TAB); 
close (CLUSTS); 
