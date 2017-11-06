#Function: Fetches the positions in the consensus where the G and A sequences were edited.
#	Steps: 	1. Receives a cluster file and extracts (a) the G and A sequence coords and (b) the edited positions in G and A sequences.
#			2. If they don't exits - creates two files per subfamily - one of G sequences and one of A sequences.
#			3. Blasts each subfamily against its consensus file (includes formatdb)
#			4. Converts the G and A positions to the position in the consensus sequence. 
#			5. Output: 	(a) Writes consensus positions to 2 separate files (G and A)
#						(b) Creates a histogram for each subfamily
#	GA_flag: 'G' or 'A' or 'GA' indicating which sequence files to fetch

#Note: 1. You need to think about what to do with L1s which are split into several consensus sequences.
# 2. additionally, in LTRs many won't have the flanking sequences, which may possibly be contained in the db sequence (needs to be checked). 

use strict;
use File::Path qw(mkpath);
use Bio::SeqIO;
( my $organism, my $class, my $family, my $pval, my $th, my $control, my $subdir, my $GA_flag ) = @ARGV;
my $dir = "../Data/$organism/$class/results/"; 
$dir .= "/$subdir" if ($subdir); 
$dir .= "/Tracks/tracks_".$organism."_".$class."_".$family."_".$pval."_".$th; 

my $repbaseDir = "../RepBase/files/RepeatDirTree_noQM"; 

#get either one (G/A) or both (G and A) seqFasta files 
opendir (DIR, $dir) || die "$dir didn't open\n";
my @files = readdir DIR; 
my $seqRegex = "seqFasta_[" . $GA_flag . "]";
my @seqFiles = grep (/$seqRegex/, @files);
my $siteRegex = "siteList_[" . $GA_flag . "]";
my @siteFiles = grep (/$siteRegex/, @files);
close(DIR);

print "@siteFiles\n"; 
print "@seqFiles\n"; 
exit; 
### Create sequence files for each subfamily (G/A specific) and get subfamily names ###
my $subfams; 
foreach my $file (@seqFiles){
	$subfams = splitFastaBySubfam($dir, $file);
	# print $file ."\n"; 
	# print "@$subfams\n";
}

### blast against the consensus sequence ### 
my @GA = split (//,$GA_flag); 
foreach my $GorA (@GA){
	my $posConsFile = $dir ."/posCons_".$GorA."_".$organism."_".$class."_".$family."_".$pval."_".$th.".txt";
	$repbaseDir .= "/$class/$family";  
	open (my $posCons_fh, ">$posConsFile") or die "$posConsFile didn't open\n";  
	foreach my $subfam (@$subfams){
		# print $subfam ."\n"; #***
		### run blast ###
		my $repbaseFile = $repbaseDir ."/". $subfam.".fa";
		my $editedFile = $dir ."/SubfamFiles/seqFasta_".$GorA."_".$subfam.".fa";  
		my $blastOutFile = $dir ."/SubfamFiles/seqFasta_".$GorA."_".$subfam."_blast.out";  
		# (my $queryFile, my $dbFile, my $outFile, my $blastType, my $cores, my $pval, my $S) = @ARGV; 
		system("blastn -query $editedFile -subject $repbaseFile -out $blastOutFile -strand plus -dust no -culling_limit 1"); #blast only plus strand
		### get consensus length ### 
		my $consStream = Bio::SeqIO->new( -file => $repbaseFile, -format => 'fasta' );
		my $consSeq = $consStream->next_seq; 
		my $length = $consSeq->length();
		
		# parse results 
		# my $subfamHist = getConsPosHist($blastOutFile, $subfam, $length); 
		# print output
		# print $posCons_fh $subfam ."\t". $@subfamHist ."\n"; #***check!
	}
	close($posConsFile);
}
### retrieve the respective positions in consensus ### 


#creates subdir in 'dir' for the output. 
#creates output files for each subfamily in the format of 'dir/subfamdir/seqFasta_[GA]_[subfam].fa'
sub splitFastaBySubfam{
	(my $dir, my $seqFile) = @_; 
	(my $GorA) = $seqFile =~ /seqFasta_([GA])/;
	my $subfamDir = $dir ."/SubfamFiles"; 
	mkpath $subfamDir;
	my %outStreams = ();
	my $inseq = Bio::SeqIO->new( -file => $dir ."/$seqFile",    -format => 'fasta' );
	while ( my $seq = $inseq->next_seq ){
		my $id = $seq->display_id(); 
		(my $subfam) = $id =~ /=([^=]+)$/; #get subfam
		if (not exists $outStreams{$subfam}){
			$outStreams{$subfam} = Bio::SeqIO->new( -file => ">".$subfamDir."/seqFasta_".$GorA."_".$subfam.".fa",  -format => 'fasta' );
		}
		$outStreams{$subfam}->write_seq($seq);
	}
	my @subfams = sort keys %outStreams;
	return \@subfams; 
}

sub getConsPosHist{
	# (my $blastOutFile, my $subfam, my $length) = @_; 
	# my @posHist = (0) x $length;
	
	# my $in = new Bio::SearchIO( -format => 'blast', -file => $blast_name );
	# while ( my $result = $in->next_result )
	# {
		# my $name = $result->query_name();
		# @indices = @{ $indices{$name} };
		##print $result->query_name(), ": @indices \n";
		# while ( my $hit = $result->next_hit )
		# {
			# $num_hits++;
			# $num_hsps = 0;
			# while ( my $hsp = $hit->next_hsp )
			# {
				# $num_hsps++;
				# my $strq      = $hsp->query_string;
				# my $strs      = $hsp->hit_string;
				# my $alignment = $hsp->homology_string;
				# my $q_start   = $hsp->start('query');
				# my $q_end   = $hsp->end('query');
				# my $s_start   = $hsp->start('subject');
				# my @sq        = split( //, $strq );
				# my @ss        = split( //, $strs );
				# my @al        = split( //, $alignment );
				# my $q_gaps    = 0;
				# my $s_gaps    = 0;

				# foreach $i ( 0 .. $#al )
				# {
					# $a = $al[$i];
					# if ( $sq[$i] eq '-' )
					# {
						# $q_gaps++;
						# next;
					# }
					# if ( $ss[$i] eq '-' )
					# {
						# $s_gaps++;
						# next;
					# }
					# $q2al[$q_start + $i - $q_gaps] = $i;
					# $al2s[$i] = $s_start + $i - $s_gaps;				
				# }			
				# foreach $ind (@indices)
				# {
					# if ($ind >= $q_start and $ind <= $q_end)
					# {
						##print "q_ind = $ind, al_ind = $q2al[$ind], s_ind = $al2s[$q2al[$ind]]\n";					
						# $hist[ $al2s[ $q2al[$ind] - 1 ] ]++;
					# }
				# }
			# }
		# }
	# }
}