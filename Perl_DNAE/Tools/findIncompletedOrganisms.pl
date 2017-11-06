#Function: creates shell-scripts to run completeBlasts.pl for all organisms that aren't finished 
#@ARGV: directory in which organisms are stored (typically "../../Data")
#Note: at this point it creates scripts for DNA (LINE LTR SINE) completion only. (see regex below)

use strict; 

my $orgDir = shift; 

my %incomplete = (); 
my %complete = (); 
opendir(ORGS, $orgDir) || print "$orgDir didn't open\n"; 
my @organisms = sort{lc($a) cmp lc($b)}(readdir(ORGS));
shift(@organisms) while ($organisms[0] =~ /^\./); 
closedir(ORGS);

foreach my $org (@organisms){
	my $classDir = $orgDir ."/". $org; 
	my @classes = &getFilesInDir($classDir); 
	foreach my $class (@classes){
		# next unless $class =~ /(LINE|LTR|SINE)/; 
		next unless $class =~ /(DNA)/; 
		my $classIncomplete = 0; 
		my $familyDir = $classDir ."/". $class . "/db"; 
		my @fams = &getFilesInDir($familyDir);
		foreach my $files_fam (@fams){
			next unless $files_fam =~ /files_(\S+)/;
			my $fam = $1; 
			my $nameDir = $familyDir."/".$files_fam; 
			my @seq_names = &getFilesInDir($nameDir);
			foreach my $seq_name (@seq_names){
				next if $seq_name =~ /\.n(hr|in|sd|si|sq|nd|ni|tm)$/;
				unless (-e "$classDir/$class/results/blasts/$fam/$seq_name".".gz"){
					$classIncomplete = 1 ;
					#print "$classDir/$class/results/blasts/$fam/$seq_name".".gz\n"; 
				}
			}
		}
		if ($classIncomplete){
			$incomplete{$org} = () unless exists $incomplete{$org}; 
			push (@{$incomplete{$org}}, $class); 
		}
		else {
			$complete{$org} = () unless exists $complete{$org}; 
			push (@{$complete{$org}}, $class); 
		}
	}
}
print "INCOMPLETE ORGANISMS:\n"; 
foreach my $org(sort(keys(%incomplete))){
	print "$org: @{$incomplete{$org}}\n"; 
}
print "********************************************************\n"; 
print "COMPLETE ORGANISMS:\n";
foreach my $org(sort(keys(%complete))){
	print "$org: @{$complete{$org}}\n"; 
}

open (SCRIPT, ">../script_completeBlasts.sh") || print "couldn't open script file\n"; 
foreach my $org(sort(keys(%incomplete))){
	print SCRIPT "nohup perl completeBlasts.pl $org 1 0 11 @{$incomplete{$org}}\n";
}
close (SCRIPT); 



sub getFilesInDir{
	(my $dir) = shift; 
	opendir (DIR, $dir) || print "couldn't open $dir\n"; 
	my @dirAr = sort{lc($a) cmp lc($b)}(readdir(DIR));
	shift(@dirAr) while ($dirAr[0] =~ /^\./);
	close(DIR); 
	return @dirAr; 
}
