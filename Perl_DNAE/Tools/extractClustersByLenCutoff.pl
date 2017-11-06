#This script must be run from one level down from "Data"
#change the location of the 
#creates a new folder inside results of Data/Organism/Class and creates new clusters inside it.
#GA_flag: flag for length filter: 0 = both, 1 = for G_seq, 2 = for A_seq; 
#NOTE: the cutoff is by the actual length and not the percent of cutoff!
use strict; 
use File::Path qw(mkpath); 

(my $organism, my $class, my $lenCutoff, my $GA_flag) = @ARGV;
my $raw_cluster_dir = "../Data/$organism/$class/results"; 
my $new_cluster_dir = "../Data/$organism/$class/results/Clusters_over_".$lenCutoff;
my %names = ();
my %nameToSerial = (); 

### copy only clusters between two sequences that are above length cutoff ###
# (file names will be same as in old directory)
my $coords; 
my $seqLen;
mkpath($new_cluster_dir); 
opendir(RAW, $raw_cluster_dir) || die "$raw_cluster_dir didn't open\n";
my @rawClustFiles = sort{lc($a) cmp lc($b)} grep {/^clusters_/}(readdir(RAW));
closedir(RAW);
#copy and extract from each file into new file
foreach my $clustFile (@rawClustFiles){
	open (OLD, "<$raw_cluster_dir/$clustFile") || die "$raw_cluster_dir/$clustFile didn't open!\n"; 
	open (NEW, ">$new_cluster_dir/$clustFile") || die "$new_cluster_dir/$clustFile didn't open!\n";
	while (my $line = <OLD>){
		if ($line =~ /^Found cluster:/){
			my $tooShort = 0; 
			my $cluster = $line; 
			while ($line = <OLD>){
				if ($line =~ /^G name = (\S+)/){
					if ($GA_flag == 0 || $GA_flag == 1){ #check length of G name
						my @seqInfo = split (/=/,$1); 
						$coords = $seqInfo[2];
						(my $start, my $end) = $coords =~ /:(\d+)-(\d+)[+-]$/;
						$seqLen = $end - $start; 
						$tooShort = 1 unless $seqLen >= $lenCutoff;
					}
				}
				if ($line =~ /^A name = (\S+)/){ #check length of A name
					if ($GA_flag == 0 || $GA_flag == 2){
						my @seqInfo = split (/=/,$1);
						$coords = $seqInfo[2]; 
						(my $start, my $end) = $coords =~ /:(\d+)-(\d+)[+-]$/;
						$seqLen = $end - $start; 
						$tooShort = 1 unless $seqLen >= $lenCutoff;
					}
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