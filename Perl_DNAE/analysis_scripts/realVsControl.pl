
#Notes: 
#	1. by default, parses only LINE, LTR, SINE, DNA, RC. (change $classRegex for other classes)
#	2. must be run from one level down from "Data". 
#	3. won't delete existing summary files. 
#	Remember to uncomment the file deletion.

use strict; 
my $classRegex = "LINE|LTR|SINE|DNA|RC"; 
#($classRegex) = @ARGV; 

my $orgDir = "../Data"; 
opendir(ORGS, $orgDir) || print "Org directory called $orgDir didn't open\n";
my @organisms = sort{lc($a) cmp lc($b)} grep{!/^\./}(readdir(ORGS));
closedir(ORGS);

my $outDir = "../DNA_editing_results"; 
open (SUMMARY, ">>$outDir/realVsControl_summary_".$$.".txt"); 
# open (REAL, ">>$outDir/real.txt") || die "couldn't open real.txt\n"; 
# open (CONTROL, ">>$outDir/control.txt") || die "couldn't open control.txt\n"; 
# open (UNCLASSIFIED, ">>$outDir/unclassified.txt") || die "couldn't open unclassified.txt\n"; 

foreach my $organism (@organisms){
	my $classDir = $orgDir . "/$organism"; 
	opendir(CLASSES, $classDir) || print "Class directory for $classDir didn't open\n";
	my @classes = sort{lc($a) cmp lc($b)} grep{!/^\./ && /$classRegex/} (readdir(CLASSES));
	closedir(CLASSES);	
	
#	print $organism . ":" ; 
#	print "@classes\n"; 
	
	foreach my $class (@classes){
		my $familyDir = "$classDir/$class/db"; 
		opendir(FAMILY, $familyDir) || print "Family directory for $familyDir didn't open\n";
		my @families = sort{lc($a) cmp lc($b)} grep{/^files_/}(readdir(FAMILY));
		closedir(FAMILY);
		for (my $i=0; $i<=$#families;$i++){ $families[$i] =~ s/files_//;} #remove 'files_' from file names - in order to fetch actual family names
			
	#	print "@families\n"; 
			
		my $clustDir = "$classDir/$class/results"; 
		opendir(CLUSTS, $clustDir) || print "Cluster directory for $clustDir didn't open\n";
		my @clustFiles = sort{lc($a) cmp lc($b)} grep{/^clusters_/ && !/control.txt$/}(readdir(CLUSTS));
		closedir(CLUSTS);
		
		#insert (real) clusterfiles to hashes, sorted by families: 
		my %cfByFams = (); 
		foreach my $cf (@clustFiles){
			$cf =~ /clusters_([^_]+)_([^_]+)_(\S+)_(1e-\d+)_(\d+)\.txt$/;  #the family is fetched with \S+ and not [^_] because some families contain underscores
			$cfByFams{$3} = () unless (exists $cfByFams{$3}); #init, if its first file for the family
			push(@{$cfByFams{$3}}, $cf); 
		}
		
		#compare file sizes of real vs control pairs(heuristics to see where strong or no signal is found)
		#init vars
		my %famSum = (); 
		my %realCount = (); 
		my %controlCount = (); 
		foreach my $family(@families){
			$famSum{$family}=0;
			$realCount{$family}=0;
			$controlCount{$family}=0;
		}
		#count for each pair if real or control file is larger;
		#Additionally, sum the size of all files of a family - to see if all files are empty (size == 0).
		foreach my $cf (@clustFiles){
			$cf =~ /clusters_([^_]+)_([^_]+)_(\S+)_(1e-\d+)_(\d+)(_control)?\.txt$/; 
			my $family = $3; 
			my $control_cf = $cf; 
			$control_cf =~ s/\.txt/_control.txt/; 
			my $realSize = -s "$clustDir/$cf"; 
			my $controlSize = -s "$clustDir/$control_cf"; 
			$famSum{$family} += $realSize + $controlSize; 
			if ($realSize > $controlSize){ #real file is larger
				$realCount{$family}++; 
			}
			elsif($realSize < $controlSize){ #control file is larger
				$controlCount{$family}++; 
			}
		}
		#print output for each family of this class
		foreach my $family (@families)
		{
			#print $organism ."\t".$family."\t" . scalar @{$cfByFams{$family}} . "\n"; #***
			my $val;
			if ($famSum{$family} == 0){
				$val = "NA"; 
				# foreach my $cf( @{$cfByFams{$family}} ){
					# unlink "$clustDir/$cf"; 
					# $cf =~ s/\.txt/_control.txt/; 
					# unlink "$clustDir/$cf"; 
				# }
			}
			elsif($realCount{$family} > $controlCount{$family}){ #real dominates 
				#2 for strong signal (in over 50% of file-pairs real > control);  1 for low signal (otherwise)
				$val = ($realCount{$family} > scalar(@{$cfByFams{$family}}) / 4) ? 2 : 1; 
			}
			elsif($realCount{$family} < $controlCount{$family}){ #control dominates
				#-2 for strong signal (in over 50% of file-pairs real < control);  -1 for low signal (otherwise)
				$val = ($controlCount{$family} > scalar(@{$cfByFams{$family}}) / 4) ? -2 : -1; 
			}
			else{ #equal score for real and control
				$val = 0; 
			}
			print SUMMARY $organism ."\t". $class ."\t". $family ."\t". $val ."\n"; 
		}
	}
	
	# $clusterFileFullPath =~ /clusters_([^_]+)_([^_]+)_([^_]+)_(1e-\d+)_(\d+)(_control)?\.txt$/;
}

# close (REAL); 
# close (CONTROL); 
# close (UNCLASSIFIED); 

close(SUMMARY); 