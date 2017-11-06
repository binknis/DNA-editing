#used for debugging an exception thrown from faulty BLAST output files

use strict;
use Bio::SearchIO;
use Bio::Root::Root;
use lib '..\Perl_for_sortClasses';
use ProcessAlignmentByLength;
use blastFormater; 
my $num_queries = 0;
my $num_hits;
my $num_hsps;
my %queries;

my $blastFile = shift; 
# my $animal = shift;
# my $virus = shift;
# my $class = shift;

my $in = new Bio::SearchIO(-format => 'blast', -file   => $blastFile);
#my $in = new Bio::SearchIO(-format => 'blast', -fh   => "STDIN");

my $progress = "progress.txt"; #create unique progress handle using PID (found in $$)

open(my $progress_handle,">$progress");
print $progress_handle "Starting...\n";
close($progress_handle);
#my $stats_name = "../Data/$animal/$virus/db/Stats_$class.txt";
#print "name = $stats_name\n";
# open(my $stats_handle,$stats_name);
# my $prob;
# while (<$stats_handle>)
# {
	# my $line = $_;
	# my $f; my $mm; my $len;
	# ($f, $mm, $len) = ($line =~ /(\S+) (\d+) (\d+)/);
	# if ($family eq $f)
	# {
		# if ($len == 0)
			# {exit();}
		# $prob = $mm / (2 * $len);
	# }
# }
print "before while in->next_result\n"; 

my $try = 1; 
while($try < 3){
	eval{
		while( my $result = $in->next_result )
		{
			print "inside\n"; 	
		}
		$try++; 
	} or do {
		print "caught\n"; 
		$try++;
		blastFormater::delZeroIdentityHSPs($blastFile);
		$in = new Bio::SearchIO(-format => 'blast', -file   => $blastFile);
	}
}
unlink($progress);





