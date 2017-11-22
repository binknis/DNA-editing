#Function: finds motifs of editing sites
#Method: for either G sequences or A sequences, it gets motifs of all sequences in a clusters output file and checks for motifs in the DB sequence
#		Uses the binomial probability
#		Only the most frequent nuc in a specific position is tested for statistical significance!!!
#$GA = if to check G or A 
#Notes: 1. sequence ID (deflines) must be identical in the DB and the siteList files
#		2. combines reusults per family (not subfamily specific)
use strict;
use Math::NumberCruncher; 
use Bio::SeqIO;

use FindBin;
use lib "$FindBin::Bin/..";  # use the parent directory of analysis_scripts directory from where this is run
use lib "$FindBin::Bin"; #because this file is in analysis_scripts, this includes that directory
my $perlDir = "$FindBin::Bin/..";

(my $dataDir, my $org, my $class, my $family, my $pvalue, my $th, my $pmotif, my $mm, my $GA , my $trackDir ) = @ARGV;

#set defaults
$GA = 'A' if ($GA eq '0' or $GA eq "");
$mm = "GA" if ($mm eq '0' or $mm eq "");
### get Nuc probs in family ###
my $nuc_name = "$dataDir/$org/$class/db/Nuc_$family.txt";
#print "name = $nuc_name\n";
open(NUC, $nuc_name) or die "no nucs: $nuc_name\n";
#create nuc probability hash: {fam} = arr of A C G T probs.
my %nuc_prob_all;
while (my $line = <NUC>)
{
	my $subfam;
	($subfam) = ($line =~ /(\S+) /);
	my $temp = $';
	my @arr = split(/ /,$temp);
	@arr = @arr[0..3];
	$nuc_prob_all{$subfam} = [@arr];		
}

### Output files ###
$pvalue = '1e-' . $pvalue unless $pvalue =~ /1e-/; 
my $suffix = $org ."_". $class ."_". $family ."_". $pvalue ."_". $th; 
my $dir = "$dataDir/".$org."/". $class ."/results";
$dir .= "/Tracks/tracks_". $suffix . "/". $mm;
$dir = $trackDir if $trackDir;
# print "my trackDr: ".$trackDir ."\n"; #***
open( my $motif_handle, ">" . $dir . "/motifs_".$GA."_".$pmotif."_".$suffix.".txt");
my $logo_file = $dir . "/logo_".$GA."_".$suffix.".txt"; 
open( my $logo_fh, ">" . $logo_file);
my $motifPerSeq_file = $dir . "/motifPerSeq_".$GA."_".$pmotif."_".$suffix.".txt"; 
open( my $mps_fh, ">" . $motifPerSeq_file); 

# print $dir . "/motifs.txt" if -e $dir . "/motifs.txt"; 
# print $dir . "/logo.txt" if -e $dir . "/logo.txt"; 
# exit;
### Get editing positions (relative to sequence start) ###
my %indices; #will contain indices of edited positions (name = unique identifier for sequence - the defline)
# open( my $sites_handle, "<". $dir . "/siteList_".$GA."_".$suffix.".txt" ) or die "siteList file didn't open for $org $class $family\n";
open( my $sites_handle, "<". $dir . "/siteList_".$GA."_".$suffix.".txt" ) or die "siteList file didn't open for $org $class $family: ".$dir . "/siteList_".$GA."_".$suffix.".txt"."\n";
while (my $line = <$sites_handle>)
{ 
	# print $line ; #***
	(my $name)  = ( $line =~ /^(\S+)\t/ );
	(my @where) = ( $'    =~ /(\d+)/g );
	(my $subfam) = $name =~ /\S+=(\S+)/;
	$subfam =~ s/\?//g;  
	$subfam =~ s/\//_/g;
	$indices{$subfam}{$name} = [@where];
}
close($sites_handle);
my $str;
my @s;
my @bf2 = ();
my @bf1 = ();
my @af1 = ();
my @af2 = ();
my $motif;

my %pos_mapping = ( 'bf2' => 0, 'bf1' => 1, 'af1' => 2, 'af2' => 3 );
my %pos_inv_mapping = ( 0 => 'bf2', 1 => 'bf1', 2 => 'af1', 3 => 'af2' );
my %nuc_mapping = ( 'x' => 0, 'a' => 1, 'c' => 2, 'g' => 3, 't' => 4 );
my %nuc_inv_mapping = (0 => 'x' , 1 => 'a' , 2 => 'c' , 3=> 'g' , 4 => 't');

my @motif_hist  = (0) x 625; #see coding explanation
my %motif2seqs = (); 
my @other_hist;
foreach my $i (0..3) { foreach my $j (0..3) { $other_hist[$i][$j] = 0; } } #4 x 4 matrix of zeros
my $count_sites = 1;
my $count = 0;

#For each subfamily, reads the DB file and finds clusters for every edited sequence (of course G or A; by $GA)
# foreach my $subfam (keys %indices)
# {  
	my $db_in_name = $dir . "/seqFasta_".$GA."_".$suffix.".fa";
	# unless(-e $db_in_name){
		# $db_in_name  = "../Data/" . $org . "/" . $class . "/db/files_" . $family . "/Seq_". $subfam;
	# }
	
	my $inseq = Bio::SeqIO->new( -file => $db_in_name,   -format => 'fasta' );
	while ( my $seq = $inseq->next_seq )
	{
		my $id = $seq->display_id();
		$id =~ s/^\d+=//; #Erase ID which isn't present in tabular output cluster files
		(my $subfam) = $id =~ /\S+=(\S+)/;
		$subfam =~ s/\?//g;  
		$subfam =~ s/\//_/g;
		if (not exists $indices{$subfam}{$id}) #skip sequences that weren't listed in edited-file
		{	
			next;
		}			
		$count++;
		$str = lc $seq->seq();
		@s = split(//,$str);
		#create arrays containing
		@bf2 = (); @bf1 = (); @af1 = (); @af2 = (); 
		foreach my $pos ( @{ $indices{$subfam}{$id} } )
		{
			$pos--; #the site pos is with 1-base start while arrays are 0-base (doesn't have anything to do with 1/0-base start in sequence IDs)
			push(@bf2, find_bf2(\@s, $pos));
			push(@bf1, find_bf1(\@s, $pos));
			push(@af1, find_af1(\@s, $pos));
			push(@af2, find_af2(\@s, $pos));			
		}			
		my @sites = @{ $indices{$subfam}{$id} };	
		for my $i (0..$#sites)
		{
			my $str = get_motif(\@s,$sites[$i]);			
			if ($str ne 'x')
			{
				print $logo_fh $str."\n";
			}
#			if ($bf2[$i] ne 'x' and $bf1[$i] ne 'x' and $af1[$i] ne 'x' and $af2[$i] ne 'x')			
#			{
#				#print $logo_fh ">$count_sites, i = $i, name = $id, site = ",$sites[$i]+1,"\n";
#				$count_sites++;
#				print $logo_fh "$bf2[$i]";
#				print $logo_fh "$bf1[$i]";
#				print $logo_fh "a";
#				print $logo_fh "$af1[$i]";
#				print $logo_fh "$af2[$i]";
#				print $logo_fh "\n";
#			}
		}
		my @probs = @{$nuc_prob_all{$subfam}};		
		#print "bf2: @bf2\n";		
		#print "bf1: @bf1\n";
		#print "af1: @af1\n";
		#print "af2: @af2\n";	
		### coding: Total of 625 cells in $motif_hist; every pos (bf2,bf1..) has 5 values (0-4).
		#			af2: values of 0,125,250,375,500. af1: 0,25,50,75,100. bf1: 0,5,10,15,20. bf2: 0-4.
		#as you can see, all ranges can't overlap others, thus any combination of values produces a unique final value. 
		$motif = 0;
		$motif += find_motif_prob('bf2', \@bf2, \@probs );			
		$motif += find_motif_prob('bf1', \@bf1, \@probs );		
		$motif += find_motif_prob('af1', \@af1, \@probs );		
		$motif += find_motif_prob('af2', \@af2, \@probs );		
		$motif_hist[$motif]++;
		$motif2seqs{$motif}{$id} = 0;
	}
# }

print $motif_handle "Motifs statistics, pvalue $pmotif, total number of sequences $count:\n";
my @sorted_hist = sort { $motif_hist[$b] <=> $motif_hist[$a] } (1..$#motif_hist); #create array of indices sorted by the values of motif_hist array
foreach my $i ( @sorted_hist )
{
	if ($motif_hist[$i] > 0)
	{
		my @str = ();
		my $j = $i;
		while ($j > 0) #build motif string
		{
			push (@str, $nuc_inv_mapping{$j % 5});
			$j /= 5; 
		}
		my @full_str;				 	
		@full_str[0..1] = @str[0..1]; $full_str[2] = $GA; @full_str[3..4] = @str[2..3]; #changed to $GA instead of GA
		print $motif_handle @full_str, ": $motif_hist[$i]\n";
		#print into motif-per-seq file (when a motif is reached it prints for all sequences with the motif)
		foreach my $idWithMotif (sort keys %{$motif2seqs{$i}}){
			print $mps_fh $idWithMotif ,"\t", @full_str, "\n";
		}
	}
}
print $motif_handle "More statistics:\n";
foreach my $i (0..3)
{
	print $motif_handle "$pos_inv_mapping{$i}: ";
	foreach my $j (0..3)
	{
		if ($other_hist[$i][$j] > 0)
		{
			print $motif_handle $nuc_inv_mapping{$j+1} , "- $other_hist[$i][$j] times, ";
		}
	}
	print $motif_handle "\n";
}
close($logo_fh); 
close($mps_fh); 
#create freqs
system("perl516 $perlDir/analysis_scripts/logoToFreq.pl $logo_file 8 2");
#***insert command to get and print normalized output

#gets the sequence in which searching for motifs will take place
# Returns: 	The G and n flanking ntds from each side.
#			If N/A (the sequences flanks aren't completely present) - returns 'x'. 
sub get_motif
{
	(my $s_ref, my $p ) = @_;	
	my @s   = @$s_ref;
	my @motif;
	my $flanking = 7; 
	if ( ($p > ($flanking-1)) and ($p < (scalar(@s)-$flanking)) )
	{
		my @inds = ($p-$flanking)..($p+$flanking);
		@motif = @s[@inds];
	}
	else
	{
		return 'x';
	}
	my $str = join('',@motif);
	$str = lc $str;
	if ($str =~ /[^acgt]/)
	{
		return 'x';
	}
	else
	{		
		return $str;	
	}	
}

sub find_bf2
{
	(my $s_ref, my $p ) = @_;	
	my @s   = @$s_ref;
	my $str = $s[ $p - 2 ];
	if ( $p > 1 and $str =~ /[acgt]/ )
	{		
		return $s[ $p - 2 ];
	}
	else
	{
		return 'x';
	}
}

sub find_bf1
{
	(my $s_ref, my $p ) = @_;
	my @s   = @$s_ref;
	my $str = $s[ $p - 1 ];
	if ( $p > 0 and $str =~ /[acgt]/ )
	{
		return $s[ $p - 1 ];
	}
	else
	{
		return 'x';
	}
}

sub find_af1
{
	( my $s_ref, my $p ) = @_;
	my @s   = @$s_ref;
	my $str = $s[ $p + 1 ];
	if ( $p < ( scalar(@s) - 1 ) and $str =~ /[acgt]/ )
	{
		return $s[ $p + 1 ];
	}
	else
	{
		return 'x';
	}
}

sub find_af2
{
	( my $s_ref, my $p ) = @_;
	my @s   = @$s_ref;
	my $str = $s[ $p + 2 ];
	if ( $p < ( scalar(@s) - 2 ) and $str =~ /[acgt]/ )
	{
		return $s[ $p + 2 ];
	}
	else
	{
		return 'x';
	}
}

sub find_motif_prob
{
	( my $where, my $arr_ref, my $nuc_prob_ref ) = @_;
	my @arr = @$arr_ref;
	my @nuc_prob = @$nuc_prob_ref;	
	my %map = ('a'=>0,'c'=>1,'g'=>2,'t'=>3);
	my %m_hist;
	$m_hist{$_}++ for (@arr);
	my $pop = ( sort { $m_hist{$b} <=> $m_hist{$a} } ( "a", "c", "g", "t" ) )[0];	
	my $n = scalar(@arr);
	my $k = $m_hist{$pop};
	my $p = $nuc_prob[$map{$pop}];	
	my $mprob = Math::NumberCruncher::Binomial($n,$k, $p); 
	#print "Here, n = $n, k = $k, p = $p, mprob = $mprob\n";
	my $what  = $pop;
	#print "where = $where, pop = $pop, $k times out of $n, prob = $mprob\n";		
	if ($mprob > $pmotif) #insignificant
	{
		$what = 'x';
	}
	else #significant
	{
		$other_hist[$pos_mapping{$where}][$nuc_mapping{$what}-1]++;
	}
	#bf1: 1-4
	my $motif = 5**$pos_mapping{$where} * $nuc_mapping{$what}; #if the prob was lower than the pmotif than this will be 0, otherwise it will move the 
	return $motif;	
}
