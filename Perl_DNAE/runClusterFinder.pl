 #Function: the file runs BLAST, finds clusters in the BLAST output and parses the cluster, for each Name (subfamily) of sequences of the organism. 
 #Input: 
 #	$allMMs = flag: 0 - check only G>A and C>T. 	1- create output for all mm combinations. 
 #If classes are specified (after all args) then retain only those classes. Otherwise, run for all classes in organisms. 
 
 use strict; 
 # use lib "/home/alu/binknis/Perl_DNAE"; 
 # use lib "$ENV{HOME}/Perl_DNAE"; 
 use FindBin;
 use lib "$FindBin::Bin";  # use the parent directory of analysis_scripts directory from where this is run
 use lib "$FindBin::Bin/analysis_scripts"; #because this file is in analysis_scripts, this includes that directory
 use removeIndexFiles; 
 use AnalyzeBlastByLength; 
 
 use Time::HiRes;
 use File::Path qw(mkpath); 
 use File::Copy;
 use IO::Compress::Gzip qw(gzip);
 use IO::Uncompress::Gunzip qw(gunzip);
 use Getopt::Long;
 use List::MoreUtils qw(uniq);
 
#Command line parameter default values
my $dataDir = $ENV{HOME} ."/Data";
my $organism = '';
my $pval_h = 2;
my $pval_l = 5;
my $th_l = 4;
my $th_h = 5;
my $cores = 0; #if set, overrides num cores for ALL parallel parameters 
my $allMMs = 1;
my $directional = 0; #analyze all 12 types of mms (e.g. GA+AG separately and not only GA) # this depends on the blast of interest
my $makeblastdb_path = ''; 
my $blastn_path = ''; 
my $useLegacyBlast = 0; 
my $more_blast_cmdline_args = ''; #e.g.: "-num_alignments 100" to restrict number of alignments (can save time)
my $blastEvalue = "1e-20"; 
my $useXMLforBlast = 0; #I read that XML BLAST output is more stable with different BioPerl releases than the text format. However, files are ~40% larger and parsing slower, hence disabled by default. 
my $override_blasts = 0; #override existing blast outfiles
my $num_parallel_blasts = 0; #set number of BLAST processes to run in parallel (see: $num_threads_blastn)
my $num_threads_blastn = 12; #number of threads per BLAST process
my $num_parallel_analyze_blast = 0; #set with number of cores if want parallel
my @classes = (); 
 
 GetOptions ("datadir|dataDir=s"  => \$dataDir,
	"organism|org=s" => \$organism,
	"pval_h=s" => \$pval_h,
	"pval_l=s" => \$pval_l,
	"th_l=i" => \$th_l,
	"th_h=i" => \$th_h,
	"allmms!" => \$allMMs,
	"directional!" => \$directional,
	
	"makeblastdb_path=s" => \$makeblastdb_path,
	"blastn_path=s" => \$blastn_path,
	"useLegacyBlast!" => \$useLegacyBlast,
	"blastEvalue=s" => \$blastEvalue,
	"usexml!" => \$useXMLforBlast,
	"blastargs=s" => \$more_blast_cmdline_args,
	
	"override_blasts" => \$override_blasts,
	
	"cores=i" => \$cores,
	"num_parallel_blasts=i" => \$num_parallel_blasts,
	"blastn_threads=i" => \$num_threads_blastn,
	"num_parallel_analyze_blast=i" => \$num_parallel_analyze_blast,
	
	"classes=s" => \@classes)
or die("Error in command line arguments\n");
 
@classes = sort(split(/,/,join(',',@classes))); #allow comma-separated list
 
 my %args = ();
$args{"dataDir"} = $dataDir; $args{"pval_h"} = $pval_h; $args{"pval_l"} = $pval_l; $args{"th_l"} = $th_l; $args{"th_h"} = $th_h; $args{"allmms"} = $allMMs;
$args{"num_parallel_analyze_blast"} = $num_parallel_analyze_blast; $args{"directional"} = $directional; 
 
$args{"bioperl_blast_read_format"} = $useXMLforBlast ? "blastxml" : "blast"; #set blast parser arg for Bioperl

if($cores > 0){ #if cores was set, need to override all parallel args (and do not enable multiple BLAST processes)
	if($num_parallel_blasts > 0){
		print "Not spawning multiple BLAST commands when 'cores' is set - multi-threading single BLAST instance\n"; 
		$num_parallel_blasts = 0; 
	}
	$num_parallel_analyze_blast = $cores;
	$num_threads_blastn - $cores; 
}
 
 
 ######################
 #######  MAIN   ######
 ######################

my $pvalue;
my $th; 
my $elapsedTime; 
my @start; 
my @end; 

my @mmDirs; 
unless($directional){
	@mmDirs = ("GA", "CT", "GC", "GT", "CA", "TA"); #GA is equivalent of AG in this algorithm
} else {
	@mmDirs = ("GA", "CT", "GC", "GT", "CA", "TA", "AG", "TC", "CG", "TG", "AC", "AT"); #GA is equivalent of AG in this algorithm
}

#create progress_$organism file
mkpath "Progress";
&write_progress($organism, "Starting ClusterFinder for $organism\n", 1); 

#Read Class names from file 
my $classDir = $dataDir. "/" . $organism; 
opendir(CLASSES, $classDir) || print "classdir didn't open\n";
my @classList = sort{lc($a) cmp lc($b)}(readdir(CLASSES));
shift(@classList) while ($classList[0] =~ /^\./); #erase "." and ".." links
closedir(CLASSES);


#Retain only classes that are listed in the classList
if (@classes){
	&getClasses(\@classList, $organism, \@classes); 
	if (scalar(@classList) == 0) { #No classes were retained - some kind of problem
		print "No classes retained in class-list for $organism in runClusterFinder.pl\n"; 
	}
}



#erase all index files where clusters will be searched for
removeIndexFiles::remove($organism, \@classList); 

###### Get list of taxa for parallelizing run #######
### Also sets up output dirs and creates cluster output files (so they exist even if none are found for the parameters)
my @iterateClass = (); 
my @iterateFam = (); 
my @iterateName = (); 

foreach my $class (@classList) {
	my $resDir = "$dataDir/$organism/$class/results"; 

	# Delete all contents of the stats file
	my $out_name = ">$dataDir/" . $organism . "/" . $class . "/db/Stats_$class.txt";
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
		
		# Create (and truncate all contents of) the cluster files (for each parameter pair)
		for(my $pval=$pval_h; $pval<=$pval_l; $pval++){
			for (my $th=$th_l; $th<=$th_h; $th++){
				my $pvalue = "1e-".$pval; 
				my $cluster_name; 
				
				#mismatch-specific dirs
				foreach my $mm (@mmDirs){
					$cluster_name = ">" . $resDir . "/" . $mm . "/clusters_" . $organism . "_" . $class . "_" . $family . "_" . $pvalue . "_" . $th . ".tab";
					open(my $c_handle, $cluster_name);
					close($c_handle);
				}
			}
		}
		
		#create family's blast output files
		mkpath "$dataDir/$organism/$class/results/blasts/$family"; 
		
		###### LOOP 3: each name in family ######
		opendir(NAME, $nameDir);
		my @nameList =  sort{lc($a) cmp lc($b)}(readdir(NAME));
		shift(@nameList) while ($nameList[0] =~ /^\./); 
		closedir(NAME);
		@nameList = grep { $_ !~ /\.n(hr|in|sd|si|sq|nd|ni|tm)$/ } @nameList; #remove index files
		
		#Save taxa iterators for runs per class-fam-subfam
		foreach my $name (@nameList){
			push(@iterateClass, $class); 
			push(@iterateFam, $family); 
			push(@iterateName, $name); 
		}
	}
}


###### BLAST all class-fam-name files to self ######
###### LOOP through all class-fam-names ######
&write_progress($organism, "Starting BLAST loops for all classes: " . join(',', uniq(@iterateClass)) . ".\n"); 
if($num_parallel_blasts > 0){ #spawn multiple BLAST processes in parallel
	my %kids; #For parallel
	my $i = -1; #For parallel
	{
	  while ($i < $#iterateName and keys %kids < $num_parallel_blasts) { #still names to run and didn't reach parallel limit (note: the class or family iterators could've also been used - same len)
		$i++; 
		my ($class, $family, $name) = ($iterateClass[$i], $iterateFam[$i], $iterateName[$i]);
		my %taxa = ("org" => $organism, "class" => $class, "fam" => $family, "name" => $name); 
		my $resDir = "$dataDir/$organism/$class/results"; 
		my $blast_archive = $resDir ."/blasts/$family/$name" . ".gz"; 
		$kids{forkBlastPerName($organism, $class, $family, $name, $dataDir, $useLegacyBlast, $makeblastdb_path, $blastn_path, $blastEvalue, $num_threads_blastn, $useXMLforBlast, $more_blast_cmdline_args, $blast_archive)} = "active";
	  }
	  {
		my $pid = waitpid(-1, 0); #wait for any child process
		if ($pid == -1) { #no child process
		  %kids = ();
		} else { #a child finished (any child)
		  delete $kids{$pid};
		}
	  }
	  redo if $i < $#iterateName or %kids;
	}
} else { #Only one BLAST process at a time
	for (my $i=0; $i<=$#iterateClass; $i++) {
		my ($class, $family, $name) = ($iterateClass[$i], $iterateFam[$i], $iterateName[$i]);
		my $resDir = "$dataDir/$organism/$class/results"; 
		
		### BLAST ###
		my $blast_archive = $resDir ."/blasts/$family/$name" . ".gz"; 
		unless (-e $blast_archive and not $override_blasts){ #BLAST file exists - uncompress it
			runBlast($organism, $class, $family, $name, $dataDir, $useLegacyBlast, $makeblastdb_path, $blastn_path, $blastEvalue, $num_threads_blastn, $useXMLforBlast, $more_blast_cmdline_args, $blast_archive);
		}			
	}		
}
&write_progress($organism, "Ended BLAST loops for all classes run: " . join(',', uniq(@iterateClass)) . ".\n"); 


###### Analyze BLAST for each name file (across all class-family-names)######
###### LOOP through all class-fam-names ######

&write_progress($organism, "Starting Analyze BLASTs for all classes: " . join(',', uniq(@iterateClass)) . ".\n"); 
if($num_parallel_analyze_blast){
	my %kids; #For parallel
	my $i = -1; #For parallel
	{
	  while ($i < $#iterateName and keys %kids < $num_parallel_analyze_blast) { #still names to run and didn't reach parallel limit (note: the class or family iterators could've also been used - same len)
		$i++; 
		my ($class, $family, $name) = ($iterateClass[$i], $iterateFam[$i], $iterateName[$i]);
		my %taxa = ("org" => $organism, "class" => $class, "fam" => $family, "name" => $name); 
		my $resDir = "$dataDir/$organism/$class/results"; 
		$kids{forkAnalyzeBlastPerName(\%taxa, \%args, $resDir)} = "active";
	  }
	  {
		my $pid = waitpid(-1, 0); #wait for any child process
		if ($pid == -1) { #no child process
		  %kids = ();
		} else { #a child finished (any child)
		  delete $kids{$pid};
		}
	  }
	  redo if $i < $#iterateName or %kids;
	}
} else { #original code - nonparallel (orginally was in the loop with BLAST)
	for (my $i=0; $i<=$#iterateName; $i++) {
		my ($class, $family, $name) = ($iterateClass[$i], $iterateFam[$i], $iterateName[$i]);
	
		### Analyze Blast with all parameters ###
		@start = Time::HiRes::gettimeofday(); #keep start time
		&write_progress($organism, "Analyzing $class $family $name ... ");
		
		my %taxa = ("org" => $organism, "class" => $class, "fam" => $family, "name" => $name); 
		
		my $formatted = AnalyzeBlastByLength::AnalyzeBlast(\%taxa, \%args); 

		if($formatted){ #save log if formatted
			print "formatted $organism $class $family $name\n" if $formatted;
		}
		
		#write analyze end time
		@end = Time::HiRes::gettimeofday(); 
		$elapsedTime = &timeDiff(\@start, \@end); 
		&write_progress($organism, "$class $family $name - Analyzed in ".$elapsedTime."\n"); 
	}
}

#erase all index files created during this run
removeIndexFiles::remove($organism, \@classList); 

######################
#######  SUBS   ######
######################
 
sub runBlast{
	my @args = @_; 
	if($#args == 0){ #array ref
		@args = @{$args[0]};
	} 
	my ($organism, $class, $family, $name, $dataDir, $useLegacyBlast, $makeblastdb_path, $blastn_path, $blastEvalue, $num_threads_blastn, $useXMLforBlast, $more_blast_cmdline_args, $blast_archive) = @args; 
	&write_progress($organism, "\t$class $family $name: Blasting ... "); 
	@start = Time::HiRes::gettimeofday(); #keep start time
	my $famDBdir = $dataDir ."/". $organism ."/". $class ."/db/files_".$family;
	
	
	my $makeblastdb_cmd; 
	my $blastn_cmd; 
	
	if($useLegacyBlast){ #Use blastall
		#set default paths for legacy blast, unless specified full-path
		$makeblastdb_path = "formatdb" unless $makeblastdb_path; 
		$blastn_path = "blastall" unless $blastn_path; 
		
		$makeblastdb_cmd = "$makeblastdb_path -i $famDBdir/$name -p F -o T"; 
		$blastn_cmd = "$blastn_path -p blastn -d $famDBdir/$name -i $famDBdir/$name -e $blastEvalue -S 1 -F F -v 0 -a $num_threads_blastn | gzip > $blast_archive";
		# print "1a. makeblastdb_cmd: $makeblastdb_cmd\n"; #***
		# print "1a. blastn_cmd: $blastn_cmd\n"; #***
	} else { #Default: Use BLAST+
		$makeblastdb_path = "makeblastdb" unless $makeblastdb_path; 
		$blastn_path = "blastn" unless $blastn_path; 
	
		$makeblastdb_cmd = "$makeblastdb_path -in $famDBdir/$name -dbtype nucl"; 
		if($useXMLforBlast){ #use XML output (else: default of text output)
			$more_blast_cmdline_args .= " -outfmt 5"; 
		}
		$blastn_cmd = "$blastn_path -query $famDBdir/$name  -db $famDBdir/$name  -evalue $blastEvalue ". ($directional ? "-strand plus" : "") . " -num_descriptions 0 -dust no -soft_masking false -num_threads $num_threads_blastn" . " " . $more_blast_cmdline_args . " | gzip > ". $blast_archive;
		# print "1b. makeblastdb_cmd: $makeblastdb_cmd\n"; #***
		# print "1b. blastn_cmd: $blastn_cmd\n"; #***
	} 
	
	#Run BLAST commands
	
	$makeblastdb_cmd =~ s/\(/\\\(/g; $makeblastdb_cmd =~ s/\)/\\\)/g; #add backslashes before "(" and ")" symbols. 
	# print "2. makeblastdb_cmd: $makeblastdb_cmd\n"; #***
	system($makeblastdb_cmd); 
	$blastn_cmd =~ s/\(/\\\(/g; $blastn_cmd =~ s/\)/\\\)/g; #add backslashes before "(" and ")" symbols. 
	# print "2. blastn_cmd: $blastn_cmd\n"; #***
	system($blastn_cmd);
	
	#write elapsed blast time
	@end = Time::HiRes::gettimeofday(); 
	$elapsedTime = &timeDiff(\@start, \@end); 
	&write_progress($organism, "$class $family $name blasted in ".$elapsedTime."\n");
}


sub forkBlastPerName{
	my $pid = fork;
	return $pid if $pid; ### parent - return PID ###
	#### child ###
	unless (defined $pid) { #fork failed
		warn "cannot fork: $!";
		return 0;
	}
	#What to do in child - START
	
	### BLAST ####
	runBlast(\@_);
	
	#What to do in child - END
	exit 0;
} 
 

sub forkAnalyzeBlastPerName{
	(my $taxa_r, my $args_r, my $resDir) = @_; 
	
	my $pid = fork;
	return $pid if $pid; ### parent - return PID ###
	#### child ###
	unless (defined $pid) { #fork failed
		warn "cannot fork: $!";
		return 0;
	}
	#What to do in child - START
	
	### Analyze Blast with all parameters ###
	# my @start = Time::HiRes::gettimeofday(); #keep start time
	# &write_progress($taxa_r->{"org"}, "Analyzing ".$taxa_r->{"name"}." ... ");
	
	my $formatted = AnalyzeBlastByLength::AnalyzeBlast($taxa_r, $args_r); 

	if($formatted){ #save log if formatted
		print "formatted ".$taxa_r->{"org"}." ".$taxa_r->{"class"}." ".$taxa_r->{"fam"}." ".$taxa_r->{"name"}."\n" if $formatted;
	}
	
	#write analyze end time
	# my @end = Time::HiRes::gettimeofday(); 
	# my $elapsedTime = &timeDiff(\@start, \@end); 
	# &write_progress($taxa_r->{"org"}, "Analyzed in ".$elapsedTime."\n"); 
	
	#What to do in child - END
	exit 0;
} 
 

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
