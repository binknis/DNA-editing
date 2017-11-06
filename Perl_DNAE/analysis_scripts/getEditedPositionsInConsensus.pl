#Function: Fetches the positions in the consensus where the G and A sequences were edited.
#	Steps: 	1. Receives a cluster file and extracts (a) the G and A sequence coords and (b) the edited positions in G and A sequences.
#			2. If they don't exits - creates two files per subfamily - one of G sequences and one of A sequences.
#			3. Blasts each subfamily against its consensus file (includes formatdb)
#			4. Converts the G and A positions to the position in the consensus sequence. 
#			5. Output: 	(a) Writes consensus positions to 2 separate files (G and A)
#						(b) Creates a histogram for each subfamily
#	$GA: 'G' or 'A' indicating which sequence files to fetch

#Note: 	1. You need to think about what to do with L1s which are split into several consensus sequences.
# 		2. additionally, in LTRs many won't have the flanking sequences, which may possibly be contained in the db sequence (needs to be checked). 
#		3. In the histogram output the first position in consensus is 1 (not 0). 
use strict;
use File::Path qw(mkpath);
use Bio::SeqIO;
use Bio::SearchIO; 
use getAll; 
( my $org, my $class, my $family, my $pval, my $th, my $control, my $subdir, my $GA ) = @ARGV;
$pval = "1e-" . $pval unless $pval =~ /1e-/; 
my $dir = "../Data/$org/$class/results"; 
$dir .= "/$subdir" if ($subdir); 
$dir .= "/Tracks/tracks_".$org."_".$class."_".$family."_".$pval."_".$th . (($control) ? "_control" : ""); 

my $subfamDir = $dir ."/SubfamFiles"; 
my $repbaseDir = "../RepBase/files/RepeatDirTree_noQM"; 
my $suffix = $org."_".$class."_".$family."_".$pval."_".$th.(($control) ? "_control" : ""); 

### Create sequence files for each subfamily (G/A specific) and get subfamily names ###
my $seqFile = "seqFasta_".$GA."_".$suffix.".fa";
my $subfams = splitFastaBySubfam($dir, $seqFile, $GA);

### get edited sites from file ###
my $siteFile = $dir."/siteList_".$GA."_".$suffix.".txt";
my $sites_ref = getEditedSites($siteFile);

### blast against the consensus sequence and create histogram of positions in consensus ### 
my $posConsFile = $dir ."/posCons_".$GA."_".$suffix.".txt";
$repbaseDir .= "/$class/$family";  
open (POSCONS, ">$posConsFile") or die "$posConsFile didn't open\n";  
foreach my $subfam (@$subfams){
	### run blast ###
	my $repbaseFile = $repbaseDir ."/". $subfam.".fa";
	my $editedFile = $subfamDir ."/seqFasta_".$GA."_".$subfam.".fa";  
	my $blastOutFile = $subfamDir ."/seqFasta_".$GA."_".$subfam."_blast.out";  
	system("perl516 Tools/runBlast.pl $editedFile $repbaseFile $blastOutFile blastn 4 20 1 1 "); #cores=4, pval=1e-20, -S=1 (only + strand), -K=1 (only 1 best hit from each region)
	# system("blastn -query $editedFile -subject $repbaseFile -out $blastOutFile -strand plus -dust no -culling_limit 1"); #blast only plus strand
	
	### parse and print results ###
	# get consensus length (needed for parsing)
	my $consStream = Bio::SeqIO->new( -file => $repbaseFile, -format => 'fasta' );
	my $consSeq = $consStream->next_seq; 
	my $consLen = $consSeq->length();
	# get and print pos hist to family
	my $subfamHist = getConsPosHist($blastOutFile, $subfam, $consLen, $sites_ref->{$subfam}); 
	print POSCONS $subfam ."\t", join(' ',@$subfamHist) ,"\n";
	# print results to specific subfamily file
	my $histFile = $subfamDir."/consPosHist_".$GA."_".$subfam.".txt";
	writeSubfamHist($histFile, $subfamHist, $consSeq);
}
close(POSCONS);

####################   		SUBROUTINES 		########################

#creates subdir in 'dir' for the output. 
#creates output files for each subfamily in the format of 'dir/subfamdir/seqFasta_[GA]_[subfam].fa'
sub splitFastaBySubfam{
	(my $dir, my $seqFile, my $GA) = @_; 
	my $subfamDir = $dir ."/SubfamFiles"; 
	mkpath $subfamDir;
	my %outStreams = ();
	my $inseq = Bio::SeqIO->new( -file => $dir ."/$seqFile",    -format => 'fasta' );
	while ( my $seq = $inseq->next_seq ){
		my $id = $seq->display_id(); 
		(my $subfam) = $id =~ /=([^=]+)$/; #get subfam
		
		if (not exists $outStreams{$subfam}){
			$outStreams{$subfam} = Bio::SeqIO->new( -file => ">".$subfamDir."/seqFasta_".$GA."_".$subfam.".fa",  -format => 'fasta' );
		}
		$outStreams{$subfam}->write_seq($seq);
	}
	my @subfams = sort keys %outStreams;
	return \@subfams; 
}

#parses a siteList file; INPUT: full path of file
sub getEditedSites{
	(my $siteFile) = @_; 
	my %sites = (); 
	my $lines = getAll::lines($siteFile);  #*** NEEDS TO BE FULL PATH!!
	die "$siteFile invalid\n" if $lines == 0; 
	foreach my $line (@$lines){
		chomp $line; 
		(my $name)  = ( $line =~ /^(\S+)\t/ );
		(my @where) = ( $'    =~ /(\d+)/g );
		(my $subfam) = ($name =~ /=([^=]+)$/); 
		$sites{$subfam}{$name} = [@where];
	}
	# foreach my $name (%sites){
		 # while (my ($k, $v) = each %{$sites{$name}}){
			# print $k ."\t" , join(' ',@$v) ,"\n"; 
		 # }
	 # }
	# exit #***
	return \%sites; 
}

sub getConsPosHist{
	(my $blastFile, my $subfam, my $consLen, my $subfamSites_ref) = @_; 
	print "$subfam\n"; #***
	my @posHist = (0) x $consLen;
	
	my $in = new Bio::SearchIO( -format => 'blast', -file => $blastFile );
	while ( my $result = $in->next_result )
	{
		my $seqID = $result->query_name();
		print $seqID . "\n"; #***
		my @indices = @{ $subfamSites_ref->{$seqID} };
		#print $result->query_name(), ": @indices \n";
		while ( my $hit = $result->next_hit )
		{
			print "hit "; 
			my @q2al; 
			my @al2s; 
			while ( my $hsp = $hit->next_hsp )
			{
				my $strq      = $hsp->query_string;
				my $strs      = $hsp->hit_string;
				my $alignment = $hsp->homology_string;
				my $q_start   = $hsp->start('query');
				my $q_end   = $hsp->end('query');
				my $s_start   = $hsp->start('subject');
				my @sq        = split( //, $strq );
				my @ss        = split( //, $strs );
				my @al        = split( //, $alignment );
				my $al_len = $hsp->length('total'); 
				my $q_gaps    = 0;
				my $s_gaps    = 0;
				@q2al = (0) x $al_len;
				@al2s = (0) x $al_len;
				foreach my $i ( 0 .. $#al )
				{
					$a = $al[$i];
					if ( $sq[$i] eq '-' )#count query gaps
					{
						$q_gaps++;
						next;
					}
					if ( $ss[$i] eq '-' )#count subject gaps
					{
						$s_gaps++;
						next;
					}
					#map query str to align str and then align str to subject str
					$q2al[$q_start + $i - $q_gaps] = $i;
					$al2s[$i] = $s_start + $i - $s_gaps;				
				}			
				foreach my $ind (@indices)
				{	
					if ($ind >= $q_start and $ind <= $q_end) #skip indices that aren't in hsp range
					{
						if ($q2al[$ind] != 0){ #skip if aligned to gap (q2al was inited to zeros)
						#print "q_ind = $ind, al_ind = $q2al[$ind], s_ind = $al2s[$q2al[$ind]]\n";					
						$posHist[ $al2s[ $q2al[$ind] - 1 ] ]++;
						}
					}
				}
			}
		}
	}
	# print $subfam.": \n"; 
	# print "@posHist\n";
	return \@posHist; 
}

sub writeSubfamHist{
	(my $histFile ,my $subfamHist, my $consSeq) = @_; 
	open( HIST, ">".$histFile ) or die "didn't open $histFile\n";
	my @str = split( //, $consSeq->seq() );
	my %mapping = ('a'=>1,'c'=>2,'g'=>3,'t'=>4);
	foreach my $i ( 0 .. $#str )
	{
		print HIST $i + 1, "\t"; #index of pos in consensus
		#amount of sites mapped to this pos
		if ( $subfamHist->[$i+1] > 0 ) 
		{
			print HIST "$subfamHist->[$i+1]\t";
		}
		else
		{
			print HIST "0\t";
		}
		#the code of nucleotide in the pos
		if ($str[$i] =~ /[acgt]/) 
		{
			print HIST $mapping{$str[$i]} , "\n";
		}
		else
		{
			print HIST "5\n";
		}
	}
	close(HIST);
}