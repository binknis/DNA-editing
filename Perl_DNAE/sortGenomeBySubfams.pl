#sort_classes.pl
#Function: Sorts genome-wide repeats to file by Classes, families and names (Taxonomy as in RepeatMasker table) 
#also creates "results" directory for each class 
#Input: 2 files - (1) 7 column data file which describes the sequences in file #2. (2) genome-wide repeats, unsorted.
#Output: Creates a directory, subdirectory and file for each Class, Family and Name respectively. 
#		Organism(host)/Class/Family/Seq_name.fasta 
#Algorithm: 1. Uses the description file to create all needed dirs and files. 
# 			2. Sorts the sequences in the sequence file into the correct dirs and files using the description file. 
#Notes: for UNIX only. 
#NOTE: output needs to be lower-case, so make sure to use $lc=1 if input fasta is upper-case 
#NOTE: "names" with question marks will be united with the same name without a question mark. The question mark will be retained in the description line! e.g. ERVK? - ERVK
#the "$newTaxa" ("$class=$family=$name") is used because some names are identical in different classes/families. 
#NOTE: only names listed in the input file will be sorted (format: 1 name per line).  

use strict; 
use File::Path qw(mkpath);
use Bio::SeqIO;
use List::Util qw(reduce);
use Statistics::Descriptive;
use getAll; 
( my $seqDataFile, my $seqFile, my $organism, my $lc, my $name_file) = @ARGV;

system("ulimit -Sn 4096"); #set the file limit to 4096, so that it doesn't fail due to "too many files open" error

### Create "familyTree", i.e. tree describing the hierarchy of Class/Family/Name of all sequences ###
open(SEQDATA,"<$seqDataFile") || die ("couldn't open description file"); 
my %familyTree = ();
my %posToTaxa = (); #chr_start_end => Class=Family=Name
my $class;
my $family;
my $name;

### create list of subfamilies to sort ###
my $subfamList = getAll::lines($name_file); 
my %namesToSort =  ();
foreach my $sfl(@$subfamList){
	$sfl =~ s/\s//g; 
	$namesToSort{$sfl} = 1;
}

while (my $line = <SEQDATA>){
	chomp; 
	next if ($line =~ /^\s*#/); #skip remark lines
	
	#create taxonomy family tree
	my @seqData = split(/\s+/,$line); 

	#save original class, family, name for sequences' defline
	#create map from position("chrom:start-end[+-]") to taxonomy ("Class=Family=Name")
	my $pos = $seqData[0] .":". $seqData[1] ."-". $seqData[2] . $seqData[3]; 
	my $taxa =  $seqData[5]."=".$seqData[6]."=".$seqData[4]; #"$class=$family=$name"
	
	#skip names that aren't in the subfam file 
	# print $seqData[4] ."\n" unless exists $namesToSort{$seqData[4]}; 
	next unless exists $namesToSort{$seqData[4]};
	
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
	mkpath "../Data/$organism/$class/results";  #create results directory
	for $family (keys %{$familyTree{$class}}){
		#create directories, files and stream for Name
		for $name (keys %{$familyTree{$class}{$family}}){
			#create directories and sequence files 
			unless (-e "../Data/$organism/$class/db/files_$family/Seq_$name"){
				mkpath "../Data/$organism/$class/db/files_$family"; 
			}
			#create, per name: 1. sequence-writing-streams, 2. Nucleotide histograms, 3. sum of all nucleotides in the name's sequences 
			my $newTaxa = "$class=$family=$name"; 
			unless(exists $seqOut{$newTaxa}){
				my $path = ">../Data/$organism/$class/db/files_$family/Seq_$name"; 
				$seqOut{$newTaxa} = Bio::SeqIO->new( -file => $path,  -format => 'fasta' );
				$hist{$newTaxa} = [ 0, 0, 0, 0 ];
			}	
		}
	}
}

#Write families to their files
#In this process: Convert sequence description lines(replaces shai's convertDB.pl)
open(my $db_in_handle, $seqFile); #maybe unneeded
my $inseq = Bio::SeqIO->new( -file => $seqFile,    -format => 'fasta' );
my %map = ('a'=>0,'c'=>1,'g'=>2,'t'=>3, 'A'=>0,'C'=>1,'G'=>2,'T'=>3); #capitals added by Binyamin to avoid need to convert later to lc
my %outCount = (); 
my %seqLen = ();
my @lenStats = ();  

while ( my $seq = $inseq->next_seq )
{
	my $id = $seq->display_id();
	($id =~ /^([^_]+)_(\S+)_(\d+)_(\d+)_([+-])\s*$/);
	my $genome = $1; 
	my $pos = $2 .":". $3 ."-". $4 . $5;
	next unless exists $posToTaxa{$pos}; #***"next" should happen only if $classFlag == 1. 
	
	#get serial number of this sequence in it's file (used for unique identifier at head of fasta defline)
	#remove ? and / characters. save original name for defline. 
	my $newTaxa = $posToTaxa{$pos}; #get "Class=Family=Name"
	$newTaxa =~ s/\?//g;  
	$newTaxa =~ s/\//_/g;
	unless (exists $outCount{$newTaxa}){
		$outCount{$newTaxa}=0; 
	}
	$outCount{$newTaxa}++; 
		
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
	
	#change seq's id, convert to lowercase and write to name's file
	my $newid = $outCount{$newTaxa}."=".$genome."=".$pos."=".$posToTaxa{$pos};
	if ($lc){
		$seq = Bio::Seq->new( -seq => lc $seq->seq(), -id  => $newid ); #create new object, with sequences converted to lower-case
	}
	else{
		$seq->display_id($newid); 
	}
	$seqOut{$newTaxa}->write_seq($seq);
}
close($db_in_handle);

#create Family "Nuc_" file, with Nucleotide occurrencies(1 line for each Name)
#note: Strand file was deleted.
#"Len_" file = list of length's of each sequence (
for $class (keys %familyTree){
	for $family (keys %{$familyTree{$class}}){
		#create family-files
		my $nuc_stats_name = ">>../Data/$organism/$class/db/Nuc_$family.txt";
		my $seq_len_name = ">>../Data/$organism/$class/db/Len_$family.txt";
		my $seq_lenStats_name = ">>../Data/$organism/$class/db/LenStats_$family.txt";
		
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
