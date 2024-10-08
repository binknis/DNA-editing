#To do: 
#If BLAST works smoothly for large SINE fams, blastFormatter::delZeroIdentityHSPs can and should be removed (BLAST used to produce some problematic output for some, usually large, files)

use strict;
use Bio::SearchIO;
use Bio::Root::Root;
use File::Path qw(mkpath); 

use FindBin;
use lib "$FindBin::Bin";  # use the parent directory of analysis_scripts directory from where this is run

use ProcessAlignmentByLength;
use blastFormatter; 

package AnalyzeBlastByLength; 

sub AnalyzeBlast {
	my $only_first_hsp = 0; #*** Flag: 0 - analyze all HSPs, 1 - only first HSP
	my $num_queries = 0;
	my $num_hits;
	my $num_hsps;
	my %queries;
	my ($taxa_r, $args_r) = @_; 
	my ($organism, $class, $family, $name) = ($taxa_r->{"org"}, $taxa_r->{"class"}, $taxa_r->{"fam"}, $taxa_r->{"name"});
	
	
	my $blastFile = $args_r->{"dataDir"} . "/$organism/$class/results/blasts/$family/$name".".gz"; 
	
	# if(-e $blastFile){
		# print $blastFile ."\n"; 
	# } elsif (){
		# print $blastFile ."\n"; 
	# } else {
		# die "zipped and non-zipped don't exist!\n";
	# }
	open(my $pipeFromZippedBlast, "gunzip -c $blastFile |") or die "open gzipped $blastFile\n"; 
	my $in = new Bio::SearchIO(-format => $args_r->{"bioperl_blast_read_format"}, -fh   => $pipeFromZippedBlast);
	my $ppid = getppid();
	my $progress = "progress_".$organism."_".$family."_".$ppid.".txt"; #create unique progress handle using PPID
	open(my $progress_handle,">>$progress");
	print $progress_handle "Name $name: Starting...\n";
	close($progress_handle);
	my $stats_name = $args_r->{"dataDir"} . "/$organism/$class/db/Stats_$family.txt";

	my $try = 1; 
	while ($try <= 2){
		  eval{
			while( my $result = $in->next_result )
			{
				$num_queries++;
				if (($num_queries % 10) == 0)
				{
					open($progress_handle,">>$progress");
					print $progress_handle "Number of queries processed = " , $in->result_count , "\n";
					close($progress_handle);
				}
				$queries{$result->query_name()} = $num_queries;
				$num_hits = 0;
				while( my $hit = $result->next_hit )
				{
					$num_hits++;
					$num_hsps = 0;
					while( my $hsp = $hit->next_hsp )
					{
						$num_hsps++;
						#my $desc = $result->query_name() . " " . $hit->name() . " ($num_hsps)"; #***unused!
						# if seq in hit was already seen as a query, don't compare current seq. with it.
						if (not exists $queries{$hit->name()})
						{
							# Avoid alignment with a sequence to itself
							if ($hsp->percent_identity() < 100)
							{
								ProcessAlignmentByLength::process_alignment_by_length(\$hsp, $result->query_name(), $hit->name(), $taxa_r, $args_r);
							}
						}
						last if $only_first_hsp; #only check first HSP for every query-subject pair
					}
				}
			}
			last; #finished successfuly - end "try" while loop
		  } or do { #erase bad format HSPs from blast file
			if ($try == 1){
				blastFormatter::delZeroIdentityHSPs($blastFile); 
				open($pipeFromZippedBlast, "gunzip -c $blastFile |") or die "open gzipped $blastFile\n"; 
				$in = new Bio::SearchIO(-format => 'blast', -fh   => $pipeFromZippedBlast);
				#save a log of classes formatted; 
				open(FORMATTED, ">>blasts_formatted_log.txt") || print "couldn't open blasts_formatted_log.txt\n"; 
				print FORMATTED $blastFile . "\n"; 
				close(FORMATTED);
				$try = 2;
			 }
			 else{ #The BLAST file was already formatted and exception was thrown (Possible cause: line of Query with all gaps). 
				#save a log of classes formated; 
				open(FORMATTED, ">>blasts_formatted_didnt_help.txt") || print "couldn't open blasts_formatted_didnt_help.txt\n"; 
				print FORMATTED $blastFile . "\n"; 
				close(FORMATTED);
				return 1; #return signal for formatting
			 }
			
		  }; 
	}
	unlink($progress);

	#returns if file was formatted
	return ($try > 1) ? 1 : 0; 
}

1; 