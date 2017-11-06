#*** SCRIPT INCOMPLETE ***#

#Function: Gets the 
#Method: for either G sequences or A sequences, it gets frequency of each nucleotide in each position upstream and downstream to the edited position. 
#$GA = if to check G or A sequence
#$normalize (flag)= 	0 - no normalization (returns count); 
#						1 - return fraction of each nuc in each pos; 
#						2 - return percentage of each nuc in each pos. 
#						3 - return fraction of each nuc in each pos + normalize by family's general nuc frequency. 
#						4 - return percentage of each nuc in each pos + normalize by family's general nuc frequency. 
#Notes: 1. possible extension: Calculate statistical significance of each result: get Nuc frequency in family; and calculate prob by binomial.
#	

#*** SCRIPT INCOMPLETE ***#

use strict;
use lib 'analysis_scripts'; 
use analysisSubs; 
use Bio::SeqIO;
use Math::Round qw(nearest);
use getAll; 

( my $org, my $class, my $family, my $pval, my $th, my $control, my $subdir, my $GA, my $range, my $normalize, my $reverse_editing ) = @ARGV;
#set defaults
$GA = 'A' if ($GA eq '0' or $GA eq "");
$pval = "1e-" . $pval unless $pval =~ /1e-/; 
### construct input filenames ###
my $dir = "../Data/" . $org . "/" . $class . "/results"; 
$dir .= "/$subdir" unless ($subdir eq 0 or $subdir eq ""); 
$dir .= "/Tracks/tracks_" . $org . "_" . $class . "_" . $family . "_" . $pval . "_" . $th;
my $siteListFile = $dir . "/siteList_".$GA."_".$org."_".$class."_".$family."_".$pval."_".$th.".txt"; #siteList_G_Lizard_LINE_L1_1e-3_3.txt
my $seqFile = $dir . "/seqFasta_".$GA."_".$org."_".$class."_".$family."_".$pval."_".$th.".fa"; #seqFasta_A_Lizard_LINE_L1_1e-3_3.fa

### get edited sites from file ###
my $sites_ref = analysisSubs::sitesFromSiteList($siteListFile);

### get nucleotide frequencies in each position ###
my @nucs = ('a', 'c', 'g', 't'); 
my @positions = (-$range .. $range);
my %nucFreq = (); 
my %posCount = (); #counter for number of positions counted for each position (isn't necessarily equal to the number of sites because some sites are close to start or end of sequence and won't have some of the positions)
foreach my $pos (@positions){
	foreach my $nuc (@nucs){
		$nucFreq{$pos}{$nuc}=0;
	}
	$posCount{$pos}=0; 
}

### fetch nucleotide frequencies around editing sites ###
my $inseq = Bio::SeqIO->new( -file => $seqFile,    -format => 'fasta' );
while ( my $seq = $inseq->next_seq ){
	my $sequence = lc $seq->seq(); #maybe del lc
	(my $nucFreq_ref) = analysisSubs::getNucFrequencyPerPos($sequence, $sites_ref->{$seq->id}, $range, $normalize, $reverse_editing); 
	for my $pos(@positions){
		for my $nuc(@nucs){
			$nucFreq{$pos}{$nuc} += $nucFreq_ref->{$pos}{$nuc}; 
		}
	}
}

### normalization (depends on $normalize arg) ###
if ($normalize > 0){
	#get nuc frequency per family
	if ($normalize > 2){
		$famNucFreq = getAll::
	}
	#sum amount of nucs counted per pos (for normalization)
	for my $pos(@positions){
		for my $nuc(@nucs){
			$posCount{$pos} += $nucFreq{$pos}{$nuc};
		}
	}
	
	#format output (percentage and round)
	my %map = ( 'a' => 0, 'c' => 1, 'g' => 2, 't' => 3);
	foreach my $pos (@positions){
		foreach my $nuc (@nucs){
			$nucFreq{$pos}{$nuc} /= $posCount{$pos} if $posCount{$pos};
			$nucFreq{$pos}{$nuc} *= $famNucFreq->[$map[$nuc]] if $normalize > 2;
			if ($normalize == 2 or $normalize == 4){
				$nucFreq{$pos}{$nuc} *= 100; #convert to percentage (if flag is '2')
				$nucFreq{$pos}{$nuc} = nearest(.01, $nucFreq{$pos}{$nuc}); #round to max of 2 digits after the decimal number
			}
		}
	}
}

### print output ###
printFreq(\%nucFreq);


#####################
######  SUBS  #######
#####################

sub printFreq{
	(my $freq_ref) = @_;
	print "position\ta\tc\tg\tt\n";
	for my $pos (sort {$a <=> $b} keys %$freq_ref){
		print $pos ."\t"; 
		for my $nuc (sort keys %{$freq_ref->{$pos}}){
			print $freq_ref->{$pos}{$nuc} . ($nuc eq 't' ? "\n" : "\t");
		}
	}
}




