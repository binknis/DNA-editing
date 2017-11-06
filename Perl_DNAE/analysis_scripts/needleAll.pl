#Function: runs needleall in parallel for one input file 
#			- Splits query file to N temp-files
#			- runs needleall from the temp-files against the subjectfile in parallel (num_processes is specified as parameter)
#			- concatenates output to single output file

#OPTION: add gzip
use strict; 
use File::Path qw(mkpath);
use getAll; 
my $dataDir; 
my $resDir; 
my $dbDir; 
my $org, my $class, my $fam, my $sf; 
my $q_file;
my $s_file;
my $out_file;
my $needleDir; 

#Check input flag - input can be either Org,class,fam,subfam (and Seq files will be fetched) 
#										or explicit input files. 
my $in_flag = shift (@ARGV); 
my $num_processes = shift (@ARGV);
if ($in_flag eq "db"){ #get files from db
	($org, $class, $fam, $sf) = @ARGV; 
	my $homeDir = $ENV{"HOME"}; 
	$dataDir =  $homeDir ."/". "Data"; 
	$resDir = "$dataDir/$org/$class/results/"; 
	$dbDir = "$dataDir/$org/$class/db/"; 
	$needleDir = $resDir . "Needles/".$fam; 
	if (-e $resDir && -e $dbDir){ #taxa exists in db
		mkpath ($needleDir); 
		$out_file = $needleDir ."/". $sf.".needle.out";
		$q_file = $s_file = $dbDir . "files_" . $fam . "/". "Seq_" . $sf; 
	}
	else { #bad input: taxa doesn't exist
		die "result file for $org $class $fam doesn't exist\n"; 
	}
}
elsif($in_flag eq "f"){ #use explicit input files
	($q_file, $s_file, $out_file) = @ARGV; 
}
else{ #bad input
	die "bad in_flag parameter - either db for database or f for explicit files\n"; 
}

##split input file into num_processes files
#create temp output file path (in needleDir and with this PID to discriminate this processes output from another's working on the same input file)
(my $q_temp_prefix) = ($q_file =~ /([^\/]+)$/);
my $pid = $$; 
$q_temp_prefix = $needleDir ."/". $q_temp_prefix . "_" . $pid; 
my $s_temp_file = $q_temp_prefix ."_subject";  
system("perl516 Tools/splitFastaFile.pl $q_file $q_temp_prefix $num_processes");
system("sed 's/:/___/g' $s_file > $s_temp_file"); 

#Run needleAll for each temp file in parallel
my @sons; 
for my $i(1 .. $num_processes){
	my $pid = fork(); 
	if ($pid){ #parent 
		push(@sons, $pid);
	}
	else{ #children
		my $q_temp_file = $q_temp_prefix .".". $i; 
		if (-e $q_temp_file){
			system("sed -i 's/:/___/g' $q_temp_file"); 
			my $temp_out_file = $q_temp_file . ".needle.out"; 
			my $cmd = "needleall -auto Y -asequence $q_temp_file -bsequence $s_temp_file -gapopen 10.0 -gapextend 0.5 -aformat score -outfile $temp_out_file"; 
			system($cmd); 
		}
		exit; 
	}
}

#wait for children and then combine files
foreach my $n (@sons){
	my $pid = waitpid($n,0);
}
#combine output files
system("sed 's/___/:/g' $q_temp_prefix.*.needle.out > $out_file"); 

#remove temp files
system("rm $q_temp_prefix.*.needle.out $s_temp_file "); 

#system("sed 's/:/___/g' $q_file | needleall -auto Y -bsequence $s_file -gapopen 10.0 -gapextend 0.5 -filter -aformat score | sed -i 's/___/:/g' > $out_file
