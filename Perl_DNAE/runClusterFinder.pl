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

 
#Command line parameter default values
my $dataDir = $ENV{HOME} ."/Data";
my $organism = '';
my $pval_h = 2;
my $pval_l = 5;
my $th_l = 4;
my $th_h = 5;
my $cores = 12;
my $allMMs = 1;
my $makeblastdb_path = ''; 
my $blastn_path = ''; 
my $useLegacyBlast = 0; 
my $more_blast_cmdline_args = ''; #e.g.: "-num_alignments 100" to restrict number of alignments (can save time)
my $blastEvalue = "1e-50"; 
my $useXMLforBlast = 0; #I read that XML BLAST output is more stable with different BioPerl releases than the text format. However, files are ~40% larger and parsing slower, hence disabled by default. 
my $override_blasts = 0; #override existing blast outfiles
my @classes = (); 
 
 GetOptions ("datadir|dataDir=s"  => \$dataDir,
	"organism|org=s" => \$organism,
	"pval_h=s" => \$pval_h,
	"pval_l=s" => \$pval_l,
	"th_l=i" => \$th_l,
	"th_h=i" => \$th_h,
	"cores=i" => \$cores,
	"allmms!" => \$allMMs,
	
	"makeblastdb_path=s" => \$makeblastdb_path,
	"blastn_path=s" => \$blastn_path,
	"useLegacyBlast!" => \$useLegacyBlast,
	"blastEvalue=s" => \$blastEvalue,
	"usexml!" => \$useXMLforBlast,
	"blastargs=s" => \$more_blast_cmdline_args,
	
	"override_blasts" => \$override_blasts,
	
	"classes=s" => \@classes)
or die("Error in command line arguments\n");
 
@classes = sort(split(/,/,join(',',@classes))); #allow comma-separated list
 
 my %args = ();
 $args{"dataDir"} = $dataDir; $args{"pval_h"} = $pval_h; $args{"pval_l"} = $pval_l; $args{"th_l"} = $th_l; $args{"th_h"} = $th_h; $args{"allmms"} = $allMMs;

 $args{"bioperl_blast_read_format"} = $useXMLforBlast ? "blastxml" : "blast"; #set blast parser arg for Bioperl
 
 
 ######################
 #######  MAIN   ######
 ######################

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

###### LOOP 1: each class in organism ######
foreach my $class (@classList) {
	my $resDir = "$dataDir/$organism/$class/results"; 
	&write_progress($organism, "$class Class:\n"); 

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
		&write_progress($organism, "\t$family Family:\n");  
		
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
		
		
		foreach my $name (@nameList) {
			next if ($name =~ /\.n(hr|in|sd|si|sq|nd|ni|tm)$/); #skip blast index files (for robustness)
			### BLAST ###
			my $blast_fileName = $resDir ."/blasts/$family/$name"; 
			my $blast_archive = $blast_fileName . ".gz"; 
			my $blastOutExists = (-e $blast_archive); 
			if ($blastOutExists and not $override_blasts){ #BLAST file exists - uncompress it
				#gunzip $blast_archive => $blast_fileName || die "gunzip failed for $blast_archive\n"; #*** changed for reading from pipe
			}
			else{ #run BLAST
				&write_progress($organism, "\t$name: Blasting ... "); 
				@start = Time::HiRes::gettimeofday(); #keep start time
				
				my $makeblastdb_cmd; 
				my $blastn_cmd; 
				
				if($useLegacyBlast){ #Use blastall
					#set default paths for legacy blast, unless specified full-path
					$makeblastdb_path = "formatdb" unless $makeblastdb_path; 
					$blastn_path = "blastall" unless $blastn_path; 
					
					$makeblastdb_cmd = "$makeblastdb_path -i $nameDir/$name -p F -o T"; 
					$blastn_cmd = "$blastn_path -p blastn -d $nameDir/$name -i $nameDir/$name -e $blastEvalue -S 1 -F F -v 0 -a $cores | gzip > $blast_archive";
					# print "1a. makeblastdb_cmd: $makeblastdb_cmd\n"; #***
					# print "1a. blastn_cmd: $blastn_cmd\n"; #***
				} else { #Default: Use BLAST+
					$makeblastdb_path = "makeblastdb" unless $makeblastdb_path; 
					$blastn_path = "blastn" unless $blastn_path; 
				
					$makeblastdb_cmd = "$makeblastdb_path -in $nameDir/$name -dbtype nucl"; 
					if($useXMLforBlast){ #use XML output (else: default of text output)
						$more_blast_cmdline_args .= " -outfmt 5"; 
					}
					$blastn_cmd = "$blastn_path -query $nameDir/$name  -db $nameDir/$name  -evalue $blastEvalue -strand plus -num_descriptions 0 -dust no -soft_masking false -num_threads $cores" . " " . $more_blast_cmdline_args . " | gzip > ". $blast_archive;
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
				&write_progress($organism, "Blasted in ".$elapsedTime."\n");
			}
			### Analyze Blast with all parameters ###
			@start = Time::HiRes::gettimeofday(); #keep start time
			&write_progress($organism, "Analyzing $name ... ");
			
			my %taxa = ("org" => $organism, "class" => $class, "fam" => $family, "name" => $name); 
			
			my $formatted = AnalyzeBlastByLength::AnalyzeBlast(\%taxa, \%args); 
	
			if($formatted){ #save log if formatted
				print "formatted $organism $class $family $name\n" if $formatted;
			}
			#compress, rename and move Blast.txt to family's results dir
			#***Not needed anymore for new blasts, see what happens for BLAST formatting
			if (not $blastOutExists and -e $blast_fileName){  
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
