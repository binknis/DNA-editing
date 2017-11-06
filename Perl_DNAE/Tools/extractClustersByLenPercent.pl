#This script must be run from one level down from "Data"
#change the location of the 
#creates a new folder inside results of Data/Organism/Class and creates new clusters inside it.
#GA_flag: flag for length filter: 0 = both, 1 = for G_seq, 2 = for A_seq; 

#Reminder: label format: 1089=hg19=chr1:145376886-145377097-=LINE=L1=HAL1

use strict; 
use File::Path qw(mkpath); 

(my $organism, my $class, my $lenPercentCutoff, my $GA_flag, my $subDir) = @ARGV;
my $raw_cluster_dir = "../Data/$organism/$class/results"; 
$raw_cluster_dir .= "/$subDir" unless $subDir eq ""; 

my $new_cluster_dir = "$raw_cluster_dir/Clusters_".$lenPercentCutoff."_percent_length_" . $GA_flag;
my %names = ();

### get names in organism's class ###
#get families
my $familyDir = "../Data/$organism/$class/db"; 
opendir(FAMILY, $familyDir) || die "Family directory for $familyDir didn't open\n";
my @familyList = sort{lc($a) cmp lc($b)} grep{/^files_/}(readdir(FAMILY));
closedir(FAMILY);

#get all names from db tree and add to the names-hash
foreach my $family (@familyList){
	my $nameDir = $familyDir . "/$family";
	opendir(NAMES, $nameDir) || die "$nameDir didn't open\n"; 
	my @namesList = grep{!/\.n(hr|in|sd|si|sq|nd|ni|tm)$/ && !/^\./} readdir(NAMES); #get all names (not index files)
	closedir(NAMES);
	foreach my $name (@namesList){
		$name =~ s/^Seq_//;
		$names{$name} = 0; 
	}
}

### get lengths of Names in Class ###
#retrieve lengths from file made from rmsk library
#create $names{name} = length
my $lengthFile = "../RepBase/files/mapping_output/name_mapping_files/db_lengths_".$class.".txt";
unless(-f $lengthFile){$lengthFile = "../RepBase/files/mapping_output/name_mapping_files/db_lengths_all.txt";} 
open (LENGTHS, "<$lengthFile") || die "$lengthFile didn't open!\n"; 
my @lengthLines = <LENGTHS>;
close(LENGTHS);

foreach my $line (@lengthLines){
	my @lenLineData = split(/\s+/,$line); #fields: 0 is Name, $# is len of consensus
	if (exists $names{$lenLineData[0]}){
		#print "exists\n"; 
		$names{$lenLineData[0]} = $lenLineData[$#lenLineData];
	}
}

#Check names to lengths.
# foreach my $name(sort keys %names){
	# print $name ."\t". $names{$name} ."\n"; 
# }
# exit; 

### copy only clusters between two sequences that are above length cutoff ###
# (file names will be same as in old directory) 
my $name; 
my $start; my $end;
opendir(RAW, $raw_cluster_dir) || die "$raw_cluster_dir didn't open\n";
my @rawClustFiles = sort{lc($a) cmp lc($b)} grep {/^clusters_/}(readdir(RAW));
closedir(RAW);
mkpath($new_cluster_dir); 
#copy and extract from each file into new file
foreach my $clustFile (@rawClustFiles){
	open (OLD, "<$raw_cluster_dir/$clustFile") || die "$raw_cluster_dir/$clustFile didn't open!\n"; 
	open (NEW, ">$new_cluster_dir/$clustFile") || die "$new_cluster_dir/$clustFile didn't open!\n";
	while (my $line = <OLD>){
		if ($line =~ /^Found cluster:/){
			my $tooShort = 0;
			my $cluster = $line; 
			while ($line = <OLD>){
				if ($line =~ /^G name = (\S+)/ && $GA_flag =~ /G/){#check length of G name
						($start, $end, $name) = ($line =~ /^G name = \d+=[^=]+=\S+:(\d+)-(\d+)[+-]=[^=]+=[^=]+=(\S+)/); #for full label: 1089=hg19=chr1:145376886-145377097-=LINE=L1=HAL1
						$tooShort = 1 if ($end-$start < $names{$name} * $lenPercentCutoff / 100); 
					#	print "checking G: " . $name ." start: ". $start . " end: " . $end . " end-start: " . ($end-$start) . " cons len:  " . $names{$name} . " bool: " . ($end-$start < $names{$name} * $lenPercentCutoff / 100) . "\n";
				}
				if ($line =~ /^A name = (\S+)/ && $GA_flag =~ /A/){ #check length of A name
						($start, $end, $name) = ($line =~ /^A name = \d+=[^=]+=\S+:(\d+)-(\d+)[+-]=[^=]+=[^=]+=(\S+)/); #for full label: 1089=hg19=chr1:145376886-145377097-=LINE=L1=HAL1
						$tooShort = 1 if ($end-$start < $names{$name} * $lenPercentCutoff / 100);
						# print "checking A: " . $name ." start: ". $start . " end: " . $end . " end-start: " . ($end-$start) . " cons len:  " . $names{$name} . " bool: " . ($end-$start < $names{$name} * $lenPercentCutoff / 100) . "\n";
				}
				$cluster .= $line;
				last if ($line =~ /^End cluster/); 
			}
			 print NEW $cluster unless ($tooShort);
		}
	}
	close(NEW); 
	close(OLD);
}