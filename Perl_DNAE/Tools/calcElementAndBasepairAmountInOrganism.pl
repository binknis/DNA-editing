#Calculates the amount of Elements and Basepairs in an organism: 
#INPUT: 1. interval file (zipped or not)
#		2. Organism (for prefix of output file and its contents)
#		3. outDir - where the 3 output files will be created
#OUTPUT: 3 files: - 1. by Class. 2. by Family. 3. by Subfamily.
#					Each with two calculated columns - #Elements and #Basepairs.
# The values will not contain the sequences and basepairs added manually to the Data directory (such as L2 for Nika)
#Note: all question-marks will be neglected and LINE and LINE? (e.g.) will be combined to LINE (for classes, fams and sfs)
use strict; 
use File::Path qw(mkpath);
(my $interval_file, my $org, my $outDir) = @ARGV; 
$outDir = "../DNA_editing_results/db_stats/AmountElementsAndBasepairs" unless $outDir; 
mkpath $outDir; 

my %seqNum = (); 
my %bpNum = (); 
my $class, my $family, my $subfam; 
my $start, my $end, my $length; 

if ($interval_file =~ /\.gz$/){ #read from pipe
	open(INTERVAL, "gunzip -c $interval_file |") or die "couldn't unzip and read $interval_file\n"; 
}
else{ #read from unzipped file
	open(INTERVAL, "<".$interval_file) or die "couldn't open $interval_file\n"; 
}

while (my $line = <INTERVAL>){
	next if $line =~ /^#/;
	chomp $line;
	$line =~ s/\?//g; #erase question marks
	my @data = split (/\s+/, $line); 
	($start, $end, $class, $family, $subfam) = ($data[1], $data[2], $data[5], $data[6], $data[4]);
	$length = $end - $start; 
	$seqNum{$class}{$family}{$subfam}++; 
	$bpNum{$class}{$family}{$subfam} += $length; 
}
close(INTERVAL); 


my $prefix = $outDir ."/". $org ."_ElementsAndBasepairs"; 
my $classFile = $prefix . "_byClass.txt"; 
my $famFile = $prefix . "_byFamily.txt"; 
my $subfamFile = $prefix ."_bySubfamily.txt";
open(C ,">".$classFile) or die "couldn't open $classFile\n"; 
open(F ,">".$famFile) or die "couldn't open $famFile\n"; 
open(SF ,">".$subfamFile) or die "couldn't open $subfamFile\n"; 

foreach my $c (sort keys %seqNum){
	my $c_seq_num = 0; 
	my $c_bp_num = 0; 
	foreach my $f (sort keys %{$seqNum{$c}}){
		my $f_seq_num=0; 
		my $f_bp_num=0; 
		foreach my $sf(sort keys %{$seqNum{$c}{$f}}){			
			print SF $org ."\t". $c ."\t". $f ."\t". $sf ."\t" . $seqNum{$c}{$f}{$sf} ."\t". $bpNum{$c}{$f}{$sf} . "\n"; 
			$f_seq_num += $seqNum{$c}{$f}{$sf}; 
			$f_bp_num += $bpNum{$c}{$f}{$sf}; 
		}
		print F $org ."\t". $c ."\t". $f ."\t". $f_seq_num ."\t". $f_bp_num . "\n"; 
		$c_seq_num += $f_seq_num; 
		$c_bp_num += $f_bp_num; 
	}
	
	print C $org ."\t". $c ."\t". $c_seq_num ."\t". $c_bp_num . "\n"; 
}
close(C); 
close(F); 
close(SF); 