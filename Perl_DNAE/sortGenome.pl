#sort_classes.pl
#Function: Sorts genome-wide repeats to file by Classes, families and names (Taxonomy as in RepeatMasker table) 
#also creates "results" directory for each class 
#Input: 2 files - (1) 7 column data file which describes the sequences in file #2. (2) genome-wide repeats, unsorted.
#Output: Creates a directory, subdirectory and file for each Class, Family and Name respectively. 
#		Organism(host)/Class/Family/Seq_name.fasta 
#Algorithm: 1. Uses the description file to create all needed dirs and files. 
# 			2. Sorts the sequences in the sequence file into the correct dirs and files using the description file. 
#Notes: for UNIX only. 
#NOTE: format of sequence files should be LOWER-CASE! (needed for cluster parsing, not for here). Set lc flag if needed. 
#NOTE: "names" with question marks will be united with the same name without a question mark. The question mark will be retained in the description line! e.g. ERVK? - ERVK
#the "$newTaxa" ("$class=$family=$name") is used because some names are identical in different classes/families. 
#Classes sorted: 
# a. Classes found in the organisms directory won't be sorted! 
# a. Specific Classes to sort can be sent as arguments (in the past they were read from classes_to_sort.txt)
# b. If no class arguments are sent, then all non-existant classes will be sorted. 
#NOTE: lower case flag added; it can be skipped; Default: converts to lower case. 


#Late edits (Starting Nov 2017): 
#1. Changed parameter input to getOpt::Long format.

#Requires a LINUX system and bedtools (change 'bedtools' call below if not in default $PATH)

#Late notes: 
#	1. If $OVERRIDE_SORTED is specified the db folder of the specified classes (or all) will be deleted. To avoid serious mistakes it currently doesn't delete the results dir (including previous blast results). Delete them manually if necessary.

use strict; 
use File::Path qw(mkpath);
use Bio::SeqIO;
use List::Util qw(reduce);
use Statistics::Descriptive;
use Getopt::Long;
use Data::Dumper;

system("ulimit -Sn 4096"); #set the file limit to 4096, so that it doesn't fail due to "too many files open" error

my $NO_TAXA_LABEL = "All"; #string to use when no class/family/name are given
my $TEMP_UNIQUE_STRING = "==="; #string used for parsing

### Commandline defaults
my $repeatIntervalFile   = '';
my $intervalFormat = ''; #currently accepts the following formats: interval|bed|rmsk (note that bed has specific format in name field or has only 4 cols chr/start/end/strand)
my $fastaFile = '';

my $genomeFile = ''; #if specified, fetches sequences (for BED or RMSK files)
my $assembly = 'NA'; #if genome is specified, this will be used to label sequences
my $DROP_ASTERISK_MARKED_IN_RMSK = 1; #if rmsk is used to fetch sequences - decides if to drop lines ending with *.

my $dataDir = $ENV{HOME} ."/Data";
my $organism = '';
my $lc = 1;
my $useFamAsName = 0; #this has higher precedence on CLASS/FAM/NAME_CONST params
my @classes = ();
my @excludeClasses = ();
my $CLASS_CONST = '';
my $FAM_CONST = '';
my $NAME_CONST = '';
my $ALL_TAXA_CONST = ''; 
my $FASTA_HEADER_REGEX = '^([^_]+)_(\S+)_(\d+)_(\d+)_([+-])\s*'; #header needs to contain assembly,chr,start,end,strand (Note: def regex previously ended with "$")
my $OPENCLOSE_PER_SEQ = 0; #to avoid issues with opening too many file handles, this can open and close the output file for each input sequence (is slower)
my $retainTempFiles = 0; 
my $tempDir = ''; #use homedir as default temp dir
my $OVERRIDE_SORTED = 0; 

#3 run options: 
#	1. Provide: $repeatIntervalFile (BED/interval), $fastaFile.
#	2. Provide: $repeatIntervalFile as (a) BED or (b) rmsk file. Must provide genome file (for seq fetching) and assembly name (for label; NA is default).
#				Details: 1)if BED, must contain chr:start-end[+-] format at beginning of name column
#	3. Provide: $fastaFile. Used when running for consensus sequences or other sequences that do not have a genome and coordinate annotation.
#	- If $dataDir isn't provided the default output dir will be used. 
#	- Provide $intervalFormat when suffix isn't self-defining .bed or .interval (e.g. rmsk or non-standard suffix like .txt)
#Example: perl516 sortGenome.pl --interval in.interval --fasta in.fa --dataDir /path/to/outdir --org orgname --lc --classes LTR,DNA
# --consttaxa Octopus

GetOptions ("interval|bed|rmsk=s" => \$repeatIntervalFile,
		  "informat=s" => \$intervalFormat,
		  "fasta|fa=s"   => \$fastaFile,
		  
		  "genome=s" => \$genomeFile,
		  "assembly=s" => \$assembly,
		  "dropAsterisks!" => \$DROP_ASTERISK_MARKED_IN_RMSK,
		  
		  "dataDir|out=s"  => \$dataDir,
		  "organism|org=s"  => \$organism,
		  "lc|tolower!"  => \$lc,
		  "classes=s" => \@classes,
		  "usefam|useFamAsName!"  => \$useFamAsName,
		  "consttaxa=s" => \$ALL_TAXA_CONST,
		  "constclass=s"  => \$CLASS_CONST,
		  "constfam=s"  => \$FAM_CONST,
		  "constname=s"  => \$NAME_CONST, 
		  "regexfasta=s"  => \$FASTA_HEADER_REGEX, 
		  "openclose|ocps!" => \$OPENCLOSE_PER_SEQ, 
		  "excludeclass=s" => \@excludeClasses, 
		  "override!" => \$OVERRIDE_SORTED,
		  
		  "tempDir" => $tempDir,
		  "retain_temp!" => $retainTempFiles)
or die("Error in command line arguments\n");

@classes = sort(split(/,/,join(',',@classes))); #allow comma-separated list
@excludeClasses = sort(split(/,/,join(',',@excludeClasses))); #allow comma-separated list

#Check interval file format
my $tempFaFile; my $tempBedFile; #needed when sequences are extracted from RMSK output
my $pseudogenome; #needed when only fasta is provided and need to create a 'pseudogenome' for downstream compatibility 

unless($intervalFormat){
	if($repeatIntervalFile =~ /.interval$/){ #interval file input
		$intervalFormat = "interval"; 
	} elsif ($repeatIntervalFile =~ /\.(bed|BED)$/) { #BED file input
		$intervalFormat = "bed"; 
		if($genomeFile and $fastaFile eq ''){
			$tempFaFile = getSeqsFromGenome($genomeFile, $repeatIntervalFile, $tempDir);
			$FASTA_HEADER_REGEX = '(\S+:\d+-\d+[+-])'; 
			$fastaFile = $tempFaFile;
			# print "fastaFile: ".$fastaFile ." tempFaFile: $tempFaFile\n"; #***
		}
	} elsif (not $repeatIntervalFile and $fastaFile) { #Fasta ONLY input (create a pseudo-genome and BED file for consistency with downstream analyses)
		$tempFaFile = $fastaFile .".".$$.".fa"; 
		$tempBedFile = $fastaFile .".".$$.".bed"; 
		$repeatIntervalFile = $tempBedFile; 
		unless($ALL_TAXA_CONST){ #This is mandatory in this option, but don't override user-specified label
			$ALL_TAXA_CONST = $NO_TAXA_LABEL;
		}
		createPsuedoGenome($fastaFile, $tempFaFile, $tempBedFile, $assembly); 
		$intervalFormat = "bed"; 
		$fastaFile = $tempFaFile;
		$FASTA_HEADER_REGEX = '^(\S+)'.$TEMP_UNIQUE_STRING.'(\S+)_(\d+)_(\d+)_([+-])\s*';
	} else { #default for non-BED and non-INTERVAL file extensions
		$intervalFormat = "rmsk"; 
		($tempBedFile, $tempFaFile) = getRMSKseqsFromGenome($repeatIntervalFile, $assembly, $tempDir);
		$repeatIntervalFile = $tempBedFile; 
		$fastaFile = $tempFaFile;
	}
}

# print "testing: " . "$repeatIntervalFile , $intervalFormat , $tempFaFile , $tempBedFile\n"; #***


### Create "familyTree", i.e. tree describing the hierarchy of Class/Family/Name of all sequences ###
open(SEQDATA,"<$repeatIntervalFile") || die ("couldn't open description file"); 
my %familyTree = ();
my %posToTaxa = (); #chr_start_end => Class=Family=Name
my $class;
my $family;
my $name;

### create list of classes to sort ###
my %classesToSort =  ();
my %classesToExclude = (); 
my %sortedClasses = (); 
#Classes were argumented - set Hash with classes to sort
if (@classes){ 
	foreach $class (@classes){
		$classesToSort{$class} = 1; 
	}
}

if (@excludeClasses){ 
	foreach $class (@excludeClasses){
		$classesToExclude{$class} = 1; 
	}
}

#get Class names which have been sorted
my @sortedClassesList = (); 
my $classDir = "$dataDir/" . $organism; 
if (-e $classDir){
	opendir(CLASSES, $classDir) || print "$classDir didn't open\n";
	@sortedClassesList = sort{lc($a) cmp lc($b)}(readdir(CLASSES));
	shift(@sortedClassesList) while ($sortedClassesList[0] =~ /^\./); #erase "." and ".." links
	closedir(CLASSES);
	foreach $class (@sortedClassesList){
		unless($OVERRIDE_SORTED){ #skip these classes
			$sortedClasses{$class} = 1;
		} elsif($OVERRIDE_SORTED and not exists $classesToExclude{$class}) { #delete sorted contents
			unless(@classes and not exists $classesToSort{$class}){
				system("rm -r $classDir/$class/db");
			}
		}
	}
}


#Map positions to taxa and create taxa family tree
while (my $line = <SEQDATA>){
	chomp $line; 
	next if ($line =~ /^\s*#/); #skip remark lines
	
	my $chr, my $start, my $end, my $strand; my $id;
	my $class, my $family, my $name; 
	my $pos; my $taxa; 
	
	#create taxonomy family tree
	my @seqData = split(/\t/,$line); 
	
	if($intervalFormat eq 'interval'){
		($chr, $start, $end, $strand, $name, $class, $family) = @seqData;  
	} elsif ($intervalFormat eq 'ucsc_rmsk'){
		die "ucsc_rmsk is still a stub\n"; 
		
	} elsif ($intervalFormat eq 'bed' or $intervalFormat eq 'rmsk'){ # expected BED format: join(/\t/, chr, start, end, coords|class|fam|name, strand)
		my $score;
		if($#seqData == 5){ #name column was specified and tentatively contains expected format (has '|')
			($chr, $start, $end, $name, $score, $strand) = @seqData; 
			if(($name =~ /\|/)){
				my @name_parts = split('\|', $name); #split by pipe
				($class, $family, $name) = @name_parts[-3..-1]; 
				# print "name_parts: " . "@name_parts\n";
				# print $class . " ", $family. " " . $name ."\n"; 
			} else {
				$class = $family = $name = $NO_TAXA_LABEL;
			}
		} elsif($#seqData == 4){ #No info specified - use NA for class family and name
			($chr, $start, $end, $strand) = @seqData;
			$class = $family = $name = $NO_TAXA_LABEL;
		} else {
			# print "@seqData\n";
			die "incompatible BED format\n";
		}
	} else {
		die "bad interval input format\n";
	}
	
	$pos = $chr .":". $start ."-". $end . $strand; #set coordinates
	
	#save original class, family, name for sequences' defline
	#create map from position("chrom:start-end[+-]") to taxonomy ("Class=Family=Name")	
	if($useFamAsName){ # for "$class=$family=$family"
		$name = $family; 
	} elsif ($ALL_TAXA_CONST ne '') { #CLASS=FAMILY=NAME are all identical
		$class = $family = $name = $ALL_TAXA_CONST; 
	} elsif ($CLASS_CONST || $FAM_CONST || $NAME_CONST){ #Use a constant label specified in in command line
		$class =  ($CLASS_CONST ? $CLASS_CONST : $class);
		$family =  ($FAM_CONST ? $FAM_CONST : $family);
		$name =  ($NAME_CONST ? $NAME_CONST : $name); 
	} else{
		$taxa =  $class."=".$family."=".$name; #"$class=$family=$name"
	}
	$taxa =  $class."=".$family."=".$name;
	
	#skip classes that aren't in class argument list to be sorted (if a list was submitted)
	my $class_no_qm = $class; 
	$class_no_qm =~ s/\?//g;
	if (@classes){ 
		next unless exists $classesToSort{$class_no_qm};
	}
	#skip classes that have been sorted already
	next if exists $sortedClasses{$class_no_qm};
	next if exists $classesToExclude{$class_no_qm};
	next if exists $posToTaxa{$pos}; #coords already exist - don't insert to familyTree (Is needed if same coords appear twice with different taxa; which happened...). 
	$posToTaxa{$pos} = $taxa; 
	
	#modify names for familyTree
	my $newTaxa = $taxa; 
	$newTaxa =~ s/\//_/g; #replace slashes with underscores (slashes are problematic for UNIX filenames).	
	$newTaxa =~ s/\?//g; #remove question marks (avoid creating new directories and/or files because of repeatMasker's uncertainty (which results in adding '?' at end of name)
	($class, $family, $name) = split (/=/, $newTaxa); 
	
	#create family tree
	$familyTree{$class}{$family}{$name}++; #***counter isn't used. 
}
close(SEQDATA); 

 
### Create files and directories for each Name ###
my %seqOut = ();
my %hist = ();

for $class (keys %familyTree){
	mkpath "$dataDir/$organism/$class/results";  #create results directory
	for $family (keys %{$familyTree{$class}}){
		#create directories, files and stream for Name
		for $name (keys %{$familyTree{$class}{$family}}){
			#create directories and sequence files 
			unless (-e "$dataDir/$organism/$class/db/files_$family/Seq_$name"){
				mkpath "$dataDir/$organism/$class/db/files_$family"; 
			}
			#create, per name: 1. sequence-writing-streams, 2. Nucleotide histograms, 3. sum of all nucleotides in the name's sequences 
			my $newTaxa = "$class=$family=$name"; 
			unless(exists $seqOut{$newTaxa}){
				my $path = "$dataDir/$organism/$class/db/files_$family/Seq_$name"; 
				if($OPENCLOSE_PER_SEQ){ 
					$seqOut{$newTaxa} = $path; #save path (don't open fh yet)
				} else {
					$seqOut{$newTaxa} = Bio::SeqIO->new( -file => ">".$path,  -format => 'fasta' ); #save file handle
				}
				$hist{$newTaxa} = [ 0, 0, 0, 0 ];
			}	
		}
	}
}


# print "dumper: " , Dumper(\%posToTaxa); #***
# if(-e $fastaFile){
	# print "fastaFile exists: $fastaFile\n"; 
# } else {
	# die "no fastaFile: $fastaFile\n";
# }

#Write families to their files
#In this process: Convert sequence description lines(replaces shai's convertDB.pl)
open(my $db_in_handle, $fastaFile) || die "can't open $fastaFile\n"; #maybe unneeded
my $inseq = Bio::SeqIO->new( -file => $fastaFile,    -format => 'fasta' );
my %map = ('a'=>0,'c'=>1,'g'=>2,'t'=>3, 'A'=>0,'C'=>1,'G'=>2,'T'=>3); #capitals added to avoid need to convert later to lc
my %outCount = (); 
my %seqLen = ();
my @lenStats = ();  
my %outputed = (); #used to avoid writing duplicate sequences to file (happens if same coords have different annotation)

while ( my $seq = $inseq->next_seq )
{
	my $genome, my $pos; 
	my @matched = ($seq->display_id() =~ /$FASTA_HEADER_REGEX/);
	if($#matched == 4){ #contains assembly: assembly_chr_start_end_strand
		$genome = $matched[0]; 
		$pos = $matched[1] .":". $matched[2] ."-". $matched[3] . $matched[4];
	} elsif ($#matched == 3) { #doesn't contain assembly: chr_start_end_strand
		$genome = $assembly; 
		$pos = $matched[0] .":". $matched[1] ."-". $matched[2] . $matched[3];
	} elsif ($#matched == 0 && $matched[0] =~ /\S+:\d+-\d+[+-]/) { #contains only pos: chr:start-end[+-]
		$genome = $assembly; 
		$pos = $matched[0]; 
	} else  {
		print "Error! var details: length: " . $#matched . "; regex: ". $FASTA_HEADER_REGEX."; display_id: ". $seq->display_id()."\n"; 
		die "bad fasta header input format\n"; 
	}
	# print "genome: $genome , pos: $pos\n"; #***
	
	next unless exists $posToTaxa{$pos}; #"next" should happen only if specific classes are analyzed
	next if exists $outputed{$pos}; #avoid writing same sequence twice to output (needed if same coords appear twice in output file; it happened...)
	$outputed{$pos}=1; 
	#get serial number of this sequence in it's file (used for unique identifier at head of fasta defline)
	#remove ? and / characters. save original name for defline. 
	my $newTaxa = $posToTaxa{$pos}; #get "Class=Family=Name"
	$newTaxa =~ s/\?//g;  
	$newTaxa =~ s/\//_/g;
	unless (exists $outCount{$newTaxa}){
		$outCount{$newTaxa}=0; 
	}
	$outCount{$newTaxa}++; 
	#change seq's id and write to name's file
	my $newid = $outCount{$newTaxa}."=".$genome."=".$pos."=".$posToTaxa{$pos};
	if ($lc){
		$seq = Bio::Seq->new( -seq => lc $seq->seq(), -id  => $newid ); #create new object, with sequences converted to lower-case
	}
	else{
		$seq->display_id($newid); 
	}
	
	if($OPENCLOSE_PER_SEQ){ #open and close fh for each sequence
		my $seqOutSingle = Bio::SeqIO->new( -file => ">>".$seqOut{$newTaxa},  -format => 'fasta' ); #save file handle (will append)
		$seqOutSingle->write_seq($seq);
		$seqOutSingle->close();
	} else { #use pre-opened fhs
		$seqOut{$newTaxa}->write_seq($seq);
	}
	
	
	#save sequence's stats
	#my $str = lc $seq->seq;	#"lc" removed by adding Capitals to %map.
	#my @s = split(//,$str); #copy sequence to array
	my @s = split(//,$seq->seq);
	foreach my $i (0..$#s) #add sequence's bases to it's class's sum
	{
		$hist{$newTaxa}[$map{$s[$i]}]++;
	}
	
	#save array of lengths for "Len_$family" file
	push ( @{$seqLen{$newTaxa}}, $seq->length()); 
}
close($db_in_handle);

#create Family "Nuc_" file, with Nucleotide occurrencies(1 line for each Name)
#note: Strand file was deleted.
#"Len_" file = list of length's of each sequence (
for $class (keys %familyTree){
	for $family (keys %{$familyTree{$class}}){
		#create family-files
		my $nuc_stats_name = ">$dataDir/$organism/$class/db/Nuc_$family.txt";
		my $seq_len_name = ">$dataDir/$organism/$class/db/Len_$family.txt";
		my $seq_lenStats_name = ">$dataDir/$organism/$class/db/LenStats_$family.txt";
		
			open( my $nuc_stats_handle, $nuc_stats_name );  
			open( my $seq_len_handle, $seq_len_name );  
			open (my $seq_lenStats_handle, $seq_lenStats_name); 
			# print $seq_lenStats_handle "#Mean\tMedian\tIQR_mean\n"; #write header in file
		#write to family-files
		for $name (keys %{$familyTree{$class}{$family}}){ #each Name - print to both files
			my $newTaxa = "$class=$family=$name"; 
			print $nuc_stats_handle "$name ";
			#sum nucleotides of all types 
			my $sum = 0;
			foreach my $i (0..3){
				$sum += $hist{$newTaxa}[$i]; 
			}
			#print histogram data (of each nucleotide)
			if ($sum == 0){
				print $seq_len_handle $organism . "$name No_sequences\n"; 
				print $seq_lenStats_handle "$name No_sequences\n";
				print $nuc_stats_handle "0 0 0 0\n"; 
			}
			else {
				foreach my $i (0..3) 
				{
					print $nuc_stats_handle $hist{$newTaxa}[$i] / $sum , " "; 
				}
				print $nuc_stats_handle "\n";
				#print lengths of each name to family-file
				print $seq_len_handle "$name @{$seqLen{$newTaxa}}\n";
				
				#get and print lengths' stats of each name to family-file
				my $lenStats_ref = &getLenStats($seqLen{$newTaxa}); 
				print $seq_lenStats_handle "$name @{$lenStats_ref}\n";
			}
		}
		close($nuc_stats_handle);
		close($seq_len_handle);
		close($seq_lenStats_handle); 
	}	
}

#Delete temporary fasta file
if(-e $tempFaFile){
	unlink($tempFaFile) unless $retainTempFiles; 
}
if(-e $tempBedFile){
	unlink $tempBedFile unless $retainTempFiles; 
}


	
	
	

sub getLenStats {
	my $lenAr_ref = shift; 
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@{$lenAr_ref});
	my @lenStats = (); 
	
	#calculate regular mean and median
	push(@lenStats, &round($stat->mean())); 
	push(@lenStats, $stat->median()); 
	#calculate IQR mean
	push(@lenStats, &round($stat->trimmed_mean(0.25))); 
	return \@lenStats; 
}

sub round{
	my $num = shift; 
	if ($num =~ /(\d+)\.(\d)/){
		if($2 >= 5){return $1+1;}
		else {return $1;} 
	}
	return $num; 
}


sub getSeqsFromGenome{
	(my $genomeFile, my $bedFile, my $tempDir) = @_; 
	my $tempFaFile = $bedFile .$$. ".temp.fa"; #temp fa file returned for sorting
	if($tempDir){
		$tempFaFile = s/(.*)\//$tempDir/;
	}
	
	system("bedtools getfasta -s -name -fi $genomeFile -bed $bedFile | sed 's/::.*//' | fold -w 80 | awk ".'\'{if($0 !~ />/) print tolower($0); else print $0;}\''." > $tempFaFile");
	# system("head -100 $tempFaFile");
	return $tempFaFile; 
}


sub getRMSKseqsFromGenome {
	my ($rmskFile, $assembly, $tempDir) = @_;
	
	my $tempBedFile = $rmskFile .$$. ".temp.bed"; #temp bed file used here
	if($tempDir){
		$tempBedFile = s/(.*)\//$tempDir/;
	}
		
	open(RMSK, "<$rmskFile") or die "couldn't open $rmskFile\n";
	open(BED, ">$tempBedFile") or die "couldn't open $tempBedFile\n";
	while(my $l = <RMSK>){
		chomp $l;  
		next if $l =~ /^(\#|\s*(SW|score))|^\s*$/; #skip comments, headers and empty lines
		if($DROP_ASTERISK_MARKED_IN_RMSK){
			next if $l =~ /\*\s*$/; 
		}
		$l =~ s/\s+/\t/g;
		$l =~ s/^\s+//;
		
		my @fields = split (/\t/, $l); #split by tab (changed from original input)
		(my $chr, my $start, my $end, my $strand, my $subfam_raw, my $classFam) = ($fields[4], $fields[5], $fields[6], $fields[8], $fields[9], $fields[10]);
		
		#modify input
		$start--; #1-base to 0-base start
		$strand =~ s/C/-/; #strand format change
		
		my $name, my $class, my $family;
		(my $name) = $subfam_raw =~ /([^|]+)$/; #extract 'name' (subfam) #Not clear why was needed
		if($classFam =~ /\//){
			($class, $family) = split('/', $classFam); 
		} else {
			$class = $family = $classFam;
		}
		
		#print interval and BED files
		my $nameForBED = join('|', $assembly."_".$chr."_".$start."_".$end."_".$strand, $class, $family, $name);
		print  BED $chr ."\t". $start."\t". $end ."\t". $nameForBED ."\t". "." ."\t". $strand ."\n";
	}
	close (BED); 
	close(RMSK);
	
	
	#Get sequences from genome
	my $tempFaFile = getSeqsFromGenome($genomeFile, $tempBedFile); 
	
	
	return ($tempBedFile, $tempFaFile); 
}


sub createPsuedoGenome {
	(my $inFasta, my $outFasta, my $outBed, my $assembly) = @_; 
	
	#open infile
	my $inseq = Bio::SeqIO->new( -file => $inFasta,    -format => 'fasta' );
	my $outseq = Bio::SeqIO->new( -file => ">".$outFasta,  -format => 'fasta' ); 
	open(my $outbed_fh,">$outBed") || die ("Error: couldn't open BED file for writing!\n"); 
	
	#create BED file and Fasta file with chromosome names
	my $chr; my $start; my $end; my $strand; 
	while ( my $seq = $inseq->next_seq ){
		# print "outbed: " . $outBed ."\n"; 
		# my $seq = Bio::Seq->new( -seq => lc $seq->seq(), -id  => $newid ); #create new object, with sequences converted to lower-case
		my $seqid = $seq->display_id(); 
		$seqid =~ s/\s+.*//; #remove whitespaces
		$seqid =~ s/=/_/g;
		$chr = $seqid; 
		$start = 0; 
		$end = $seq->length()-1; 
		$strand = '+';
		$seqid = $assembly .$TEMP_UNIQUE_STRING. join('_', $chr, $start, $end, $strand); 
		$seq->display_id($seqid); 
		$seq->desc('');
		$outseq->write_seq($seq);
		print $outbed_fh join("\t",$chr, $start, $end, join('|', $NO_TAXA_LABEL, $NO_TAXA_LABEL, $chr), ".", "+") , "\n"; #BED: chr, start, end, name(seq_header), ., strand
	}

	#close files
	$inseq->close();
	$outseq->close();
	close($outbed_fh); 
}




