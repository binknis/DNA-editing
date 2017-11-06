 #Function: the file runs BLAST, finds clusters in the BLAST output and parses the cluster, for each Name (subfamily) of sequences of the organism. 
 #Input: 
 #	$allMMs = flag: 0 - check only G>A and C>T. 	1- create output for all mm combinations. 
 #If classes are specified (after all args) then retain only those classes. Otherwise, run for all classes in organisms. 
 
 use strict; 
 use lib "/home/alu/binknis/Perl_DNAE"; 
 use lib "$ENV{HOME}/Perl_DNAE"; 
 use Time::HiRes;
 use removeIndexFiles; 
 use File::Path qw(mkpath); 
 use File::Copy;
 use IO::Compress::Gzip qw(gzip);
 use IO::Uncompress::Gunzip qw(gunzip);
 use AnalyzeBlastByLength; 
 
 ######################
 #######  MAIN   ######
 ######################
my $num_to_splice = 7; 
(my $organism, my $pval_h, my $pval_l, my $th_l, my $th_h, my $cores, my $allMMs, my $tabular) = splice(@ARGV,0,$num_to_splice);
my @argClasses = @ARGV;
if ($allMMs ne '0' && $allMMs ne '1'){ #no $allMMs inputed, rather a class was accidently extracted from ARGV list
	push (@argClasses, $allMMs); 
	$allMMs = 1; #DEFAULT - check all mismatches
}
if ($tabular ne '0' && $tabular ne '1'){ #no $tabular inputed, rather a class was accidently extracted from ARGV list
	push (@argClasses, $tabular); 
	$tabular = 1; #DEFAULT - print in tabular output & each file in mismatch-specific subdir
}
my $pvalue;
my $th; 
my $elapsedTime; 
my @start; 
my @end; 

$cores = 1 if ($cores == 0); #set 1 as default for core amount for blasts. 
my @mmDirs = ("GA", "CT", "GC", "GT", "CA", "TA"); #GA is equivalent of AG in this algorithm

#create progress_$organism file
mkpath "Progress";
&write_progress($organism, "Starting ClusterFinder for $organism\n", 1); 

#Read Class names from file 
my $classDir = "../Data/" . $organism; 
opendir(CLASSES, $classDir) || print "classdir didn't open\n";
my @classList = sort{lc($a) cmp lc($b)}(readdir(CLASSES));
shift(@classList) while ($classList[0] =~ /^\./); #erase "." and ".." links
closedir(CLASSES);

#Retain only classes that are listed in the classFile (if classFlag == 1)
if (@argClasses){
	&getClasses(\@classList, $organism, \@argClasses); 
	if (scalar(@classList) == 0) { #No classes were retained - some kind of problem
		print "No classes retained in class-list for $organism in runClusterFinder.pl\n"; 
	}
}
#erase all index files where clusters will be searched for
removeIndexFiles::remove($organism, \@classList); 

###### LOOP 1: each class in organism ######
foreach my $class (@classList) {
	my $resDir = "../Data/$organism/$class/results"; 
	&write_progress($organism, "$class Class:\n"); 

	# Delete all contents of the stats file
	my $out_name = ">../Data/" . $organism . "/" . $class . "/db/Stats_$class.txt";
	open(my $temp,$out_name);
	close($temp);
	
	# Create result dirs
	foreach my $mm(@mmDirs){
		mkpath $resDir . "/$mm"; 
	}
	
	###### LOOP 2: each family in class ######
	my $familyDir = $classDir ."/$class/db"; 
	opendir(FAMILY, $familyDir) || print "Family directory for $familyDir didn't open\n";
	my @familyList = sort{lc($a) cmp lc($b)}(readdir(FAMILY));
	shift(@familyList) while ($familyList[0] =~ /^\./); 
	closedir(FAMILY);
	
	foreach my $family (@familyList) {
		next unless (-d $familyDir ."/$family" ); #skip non-directory files
		my $nameDir = $familyDir . "/$family";
		next if ($family !~ /files_(\S+)/); #skip directories that don't contain sequences
		$family = $1; #erase "files_" prefix from family-name
		&write_progress($organism, "\t$family Family:\n");  
		
		# Create (and truncate all contents of) the cluster files (for each parameter pair)
		for(my $pval=$pval_h; $pval<=$pval_l; $pval++){
			for (my $th=$th_l; $th<=$th_h; $th++){
				my $pvalue = "1e-".$pval; 
				my $cluster_name; 
				
				#mismatch-specific dirs
				foreach my $mm (@mmDirs){
					next if ( not $tabular and ($mm eq 'GA' or $mm eq 'CT') ); #Don't create subdirs for GA and CT output when non-tabular output format was chosen (tabular=0)
					$cluster_name = ">" . $resDir . "/" . $mm . "/clusters_" . $organism . "_" . $class . "_" . $family . "_" . $pvalue . "_" . $th . ($tabular ? ".tab" : ".txt");
					open(my $c_handle, $cluster_name);
					close($c_handle);
				}
				
				unless ($tabular){ #not tabular output - create files in "results" dir for GA and CT (old)
					#real (GA)
					$cluster_name = ">" . $resDir . "/clusters_" . $organism . "_" . $class . "_" . $family . "_" . $pvalue . "_" . $th . ".txt";
					open(my $c_handle, $cluster_name);
					close($c_handle);
					#control (CT)
					$cluster_name = ">" . $resDir . "/clusters_" . $organism . "_" . $class . "_" . $family . "_" . $pvalue . "_" . $th . "_control.txt";
					open(my $c_handle, $cluster_name);
					close($c_handle);
				}
			}
		}
		
		#create family's blast output files
		mkpath "../Data/$organism/$class/results/blasts/$family"; 
		
		###### LOOP 3: each name in family ######
		opendir(NAME, $nameDir);
		my @nameList =  sort{lc($a) cmp lc($b)}(readdir(NAME));
		shift(@nameList) while ($nameList[0] =~ /^\./); 
		closedir(NAME);
		foreach my $name (@nameList) {
			next if ($name =~ /\.n(hr|in|sd|si|sq|nd|ni|tm)$/); #skip blast index files (for robustness)
			### BLAST ###
			my $blast_fileName = $resDir ."/blasts/$family/$name"; 
			my $blast_archive = $blast_fileName . ".gz"; 
			my $doNotBlast = (-e $blast_archive); 
			if ($doNotBlast){ #BLAST file exists - uncompress it
				#gunzip $blast_archive => $blast_fileName || die "gunzip failed for $blast_archive\n"; #*** changed for reading from pipe
			}
			else{ #run BLAST
				&write_progress($organism, "\t$name: Blasting ... "); 
				@start = Time::HiRes::gettimeofday(); #keep start time
				
				my $formatdb = "formatdb -i $nameDir/$name -p F -o T"; 
				$formatdb =~ s/\(/\\\(/g; $formatdb =~ s/\)/\\\)/g; #add backslashes before "(" and ")" symbols. 
				system($formatdb); 
				my $blast = "blastall -p blastn -d $nameDir/$name -i $nameDir/$name -e 1e-50 -S 1 -F F -v 0 -a $cores > $blast_fileName";
				$blast =~ s/\(/\\\(/g; $blast =~ s/\)/\\\)/g; #add backslashes before "(" and ")" symbols. 
				system($blast);
				
				#write elapsed blast time
				@end = Time::HiRes::gettimeofday(); 
				$elapsedTime = &timeDiff(\@start, \@end); 
				&write_progress($organism, "Blasted in ".$elapsedTime."\n");
			}
			### Analyze Blast with all parameters ###
			@start = Time::HiRes::gettimeofday(); #keep start time
			&write_progress($organism, "Analyzing $name ... ");
			
			# my $analyzeBlast = "perl AnalyzeBlastByLength.pl $name $organism $class $family $pval_h $pval_l $th_l $th_h"; 
			# $analyzeBlast =~ s/\(/\\\(/g; $analyzeBlast =~ s/\)/\\\)/g; #add backslashes before "(" and ")" symbols. 
			# system($analyzeBlast);
			#the following line should replace the preceding 3
			my $formatted = AnalyzeBlastByLength::AnalyzeBlast($name, $allMMs, $organism, $class, $family, $pval_h, $pval_l, $th_l, $th_h, $tabular); 
			if($formatted){ #save log if formatted
				print "formatted $organism $class $family $name\n" if $formatted;
			}
			#compress, rename and move Blast.txt to family's results dir
			#(Needed only if BLAST was created now or it was formatted)
			if (not $doNotBlast){  
				gzip $blast_fileName => $blast_archive or print "gzip failed for $organism $class $family $blast_fileName\n";
			}
			unlink($blast_fileName);
			
			#write analyze end time
			@end = Time::HiRes::gettimeofday(); 
			$elapsedTime = &timeDiff(\@start, \@end); 
			&write_progress($organism, "Analyzed in ".$elapsedTime."\n"); 
		}
	}
}

#erase all index files created during this run
removeIndexFiles::remove($organism, \@classList); 

######################
#######  SUBS   ######
######################
 
sub getClasses{
	my $classes_ref = shift; 
	my $organism = shift; 
	my $argClassList = shift; 
	my %classesToRetain = ();  
	
	if ($argClassList){ #specific classes were argumented to command line
		foreach my $class(@{$argClassList}){
			$classesToRetain{$class} = 1; 
		}
	} 
	else{ #no argumented classes - use classes.txt class (this was the initial method, the "if" is the new feature using argumented classes in command line)
		open(my $classes_fh, "classes.txt") || print "Couldn't open classes.txt in runClusterFinder.pl - Running for all classes of $organism\n";
		while (<$classes_fh>){
			foreach my $class (split){
				$classesToRetain{$class} = 1; 
			}
		}
		close($classes_fh); 
	}
	
	#create new classList (intersection of classes in classFile and in classesToRetain)
	my @newClassList = (); 
	foreach my $class (@{$classes_ref}){
		if (exists $classesToRetain{$class}){
			push(@newClassList, $class); 
		}
	}
	@{$classes_ref} = @newClassList;
}
 

#Write the progress-log to progress_$organism file
sub write_progress{
my $organism = shift; 
my $toPrint = shift;
my $flag = shift; 
	if ($flag){
		open (PROGRESS, ">Progress/progress_$organism"."_".$$); #"$$" holds PID
	}else{
		open (PROGRESS, ">>Progress/progress_$organism"."_".$$); 
	}
	print PROGRESS $toPrint;
	close(PROGRESS); 
}

#Calculate the difference between two time points (saved as array[0] = seconds, array[1] = microseconds)
sub timeDiff{
	my $start_ref = shift; 
	my @start = @{$start_ref}; 
	my $end_ref = shift; 
	my @end = @{$end_ref}; 
	#calculate difference (seconds, microseconds)
	my @diff;
	$#diff = 2; 
	for (my $i=0; $i<2;$i++){
		@diff[$i] = @end[$i] - @start[$i]; 
	}

	#convert to Hour:Min:Sec.fraction_of_sec format
	my $hour = int($diff[0] / 3600); 
	$diff[0] -= $hour * 3600; 
	my $min = int($diff[0] / 60); 
	$diff[0] -= $min * 60; 
	my $sec = $diff[0];
	my $frac = int(($diff[1] / 1000000) * 10); #save only tenthes of the second
	if ($frac < 0) {
		$frac = 10 + $frac;
		$sec--; 
	} 

	#keep two digit format
	if ($hour<10){$hour = 0 . $hour;}
	if ($min<10){$min = 0 . $min;}
	if ($sec<10){$sec = 0 . $sec;}
	return "$hour:$min:$sec.$frac\n";
}