#Function: Compares the occurrence of a substring in edited vs. unedited sequences
#Input:
# 1. BED file with edited elements. 
# 2. Fasta file with edited elements and background 

#Output: 
#1. Creates a sequence file with all 'control' (unedited) elements. 

#Note: uses getAll.pm, therefore must run from level under "Data".
use strict; 
use File::Path qw(mkpath);
use Bio::SeqIO;
use lib $ENV{HOME} . "/Perl_DNAE";
use getAll; 
use lib $ENV{HOME} . "/Perl_DNAE" . "/analysis_scripts";
use analysisSubs; 

# use Data::Dumper; #***

(my $org, my $class, my $orgSeqFasta, my $uneditedOutfile) = @ARGV; 
$org = "Zebrafinch" unless $org; 
$class = "LTR" unless $class; 
#Used 
$orgSeqFasta = "/home/alu/binknis/Rscripts/DNAE/Choose_best_params/lowConfidence_byClusters_AmoreDiv_mostGMapGinPairs_bestSources/perOrgFiles/Zebrafinch/seqFasta_bothGandA_Zebrafinch_LTR_allFams_1e-0_5.fa" unless $orgSeqFasta;
#This file is used for calculating "edited" frequencies
my $orgSeqFasta_GseqsOnly = "/home/alu/binknis/Rscripts/DNAE/Choose_best_params/lowConfidence_byClusters_AmoreDiv_mostGMapGinPairs_bestSources/perOrgFiles/Zebrafinch/seqFasta_G_Zebrafinch_LTR_allFams_1e-0_5.lowerCase.fa";

unless($uneditedOutfile){
	$uneditedOutfile = $orgSeqFasta; #where all the unedited files will be writted to - for control
	$uneditedOutfile =~ s/.fa/.control.fa/;
}
my $nuc = "g"; #nucToAnalyze
my @nucs = ('a', 'c', 'g', 't'); 
#Get list of subfams and coords
(my $subfams, my $coords) = getSubfamsAndCoords($orgSeqFasta); 

createConrolFileFromDB($org, $class, $subfams, $coords, $uneditedOutfile); 

### Get nucleotide frequencies
#edited
my $range = 1; #get only frequencies at +1
my $retfmt = 0; #return frequencies
my $nfe =  analysisSubs::getNucFreqFromFasta($orgSeqFasta_GseqsOnly, $range, $retfmt); #get nuc frequency for edited elements
print "Edited:\n"; 
analysisSubs::printFreqHash($nfe->{$nuc}); 
my $total = 0; 
foreach my $n (@nucs){
	$total += $nfe->{$n}{"0"}{$n}; 
}
print "Total seq lengths: $total\n"; 

#control
my $nfue = analysisSubs::getNucFreqFromFasta($uneditedOutfile, $range, $retfmt); #get nuc frequency for unedited elements
print "Unedited background:\n"; 
analysisSubs::printFreqHash($nfue->{$nuc}); 
my $total = 0; 
foreach my $n (@nucs){
	$total += $nfue->{$n}{"0"}{$n}; 
}
print "Total seq lengths: $total\n"; 

#Get two hashes from a fasta file, one with subfams and another with coords
#These are used to select the correct sequences for the 'control' file
#	i.e. only sequences in the edited-subfams but not the edited seqs are selected.
sub getSubfamsAndCoords{
	(my $seqFile) = @_; 
	my %subfamHash = (); 
	my %coordHash = (); 
	my $inseq = Bio::SeqIO->new( -file => $seqFile,    -format => 'fasta' );
	while ( my $seq = $inseq->next_seq ){
		my $id = $seq->display_id(); 
		(my $coords, my $subfam) = $id =~ /=([^=]+:\d+-\d+[+-])=\S+=([^=]+)$/; #get coords from id
		$subfamHash{$subfam} = 0; 
		$coordHash{$coords} = 0;
	}
	# print Dumper(\%subfamHash); #***
	# print Dumper(\%coordHash); #***
	return (\%subfamHash, \%coordHash); 
}

sub createConrolFileFromDB{
	(my $org, my $class, my $sfs_ref, my $coords_ref, my $outfile) = @_; 
	#open stream to output file
	my $outStream = Bio::SeqIO->new( -file => ">".$outfile,  -format => 'fasta' ); 
	#get wanted sequences
	my $fams = getAll::families($org, $class); 
	foreach my $fam (@$fams){
		my $nameFiles = getAll::name_files($org, $class, $fam, 1); 
		foreach my $nf (@$nameFiles){
			(my $sf) = $nf =~ /.*\/Seq_(\S+)/; 
			# print $sf ."\n";
			if(exists $sfs_ref->{$sf}){ #subfam is in hash (contains edited elements)
				# system("cat $nf >> $outfile") #used at first to select all
				my $inseq = Bio::SeqIO->new( -file => $nf,    -format => 'fasta' );
				while ( my $seq = $inseq->next_seq ){
					my $id = $seq->display_id(); 
					(my $coords) = $id =~ /=([^=]+:\d+-\d+[+-])=/; #get coords from id
					# print $coords ."\n"; 
					$seq = Bio::Seq->new( -seq => lc $seq->seq(), -id  => $seq->display_id()); #convert sequence to lowercase
					unless (exists $coords_ref->{$coords}){ #write the sequence to output if the sequence's coords DO NOT exist in the coords list (get only unedited elements)
						$outStream->write_seq($seq);
					}
				}
			}
			else{ #subfamily isn't in hash (didn't contain edited elements)
				next; #don't get the sequences
			}
			
		}
		
	}
}
