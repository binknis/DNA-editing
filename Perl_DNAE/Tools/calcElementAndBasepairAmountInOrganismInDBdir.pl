#Calculates the amount of Elements and Basepairs in an organism: 
#INPUT: 1. Organism (for prefix of output file and its contents)
#		2. outDir - where the 3 output files will be created
#OUTPUT: 3 files: - 1. by Class. 2. by Family. 3. by Subfamily.
#					Each with two calculated columns - #Elements and #Basepairs.
use strict; 
my $home = $ENV{"HOME"}; 
use lib $home ."/Perl_DNAE"; 
use getAll; 
use Bio::SeqIO;
use File::Path qw(mkpath);
(my $org, my $outDir) = @ARGV; 

$outDir = $home ."/DNA_editing_results/db_stats/AmountElementsAndBasepairs" unless $outDir; 
mkpath $outDir; 

my %seqNum = (); 
my %bpNum = (); 
my $class, my $family, my $subfam; 
my $start, my $end, my $length; 

my $classes = getAll::classes($org); 
foreach my $class (@$classes){
	my $fams = getAll::families($org, $class); 
	foreach my $family (@$fams){
		my $subfams = getAll::names($org, $class, $family); 
		foreach my $subfam (@$subfams){
			my $subfam_file = $home ."/Data/". $org ."/". $class ."/db/files_".$family."/Seq_".$subfam;  
			(my $seqs, my $bps) = getSeqAndBpCount($subfam_file); 
			$seqNum{$class}{$family}{$subfam} = $seqs;
			$bpNum{$class}{$family}{$subfam} = $bps; 
		}
	}
}

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



#### SUBS ####
sub getSeqAndBpCount{
	(my $sf_file) = @_; 
	my $in = Bio::SeqIO->new( -file => $sf_file, -format => 'fasta' );
	my $seqs=0; 
	my $bps=0; 
	while (my $seq = $in->next_seq){
		$seqs++; 
		$bps += $seq->length; 
	}
	return ($seqs, $bps); 
}