#Function: creates shell-scripts to run completeBlasts.pl for all organisms that aren't finished 
#dataDir: directory in which organisms are stored (typically HOME ."/Data")
#Note: at this point it creates scripts for LTR completion only. (see regex below)

use strict; 

(my @classes_in) = @ARGV; 
my $dataDir = $ENV{HOME} ."/Data"; 
my %incomplete = (); 
my %complete = ();  		
opendir(ORGS, $dataDir) || print "$dataDir didn't open\n"; 
my @organisms = sort{lc($a) cmp lc($b)}(readdir(ORGS));
shift(@organisms) while ($organisms[0] =~ /^\./); 
closedir(ORGS);

my $class_regex; 
if (@classes_in){ #classes inputed - check only them
	$class_regex = "(" . join ("|", @classes_in) . ")"; 
}
else { #no class input - check all
	$class_regex = '\S+';  
}

foreach my $org (@organisms){
	my $classDir = $dataDir ."/". $org; 
	my @classes = &getFilesInDir($classDir); 
	foreach my $class (@classes){
		# next unless $class =~ /(LINE|LTR|SINE)/; 
		# next unless $class =~ /LTR/; 
		next unless $class =~ /$class_regex/; 
		my $classIncomplete = 1; 
		my $resultDir = $classDir ."/". $class . "/results"; 
		my @results = &getFilesInDir($resultDir);
		
		foreach my $result_file (@results){
			if ($result_file =~ /^clusters_/){
				$classIncomplete = 0; 
				last;
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
	print "$org @{$incomplete{$org}}\n"; 
}
print "********************************************************\n"; 
print "COMPLETE ORGANISMS:\n";
foreach my $org(sort(keys(%complete))){
	print "$org @{$complete{$org}}\n"; 
}

# open (SCRIPT, ">../script_completeBlasts.sh") || print "couldn't open script file\n"; 
# foreach my $org(sort(keys(%incomplete))){
	# print SCRIPT "nohup perl completeBlasts.pl $org 1 0 11 @{$incomplete{$org}}\n";
# }
# close (SCRIPT);



sub getFilesInDir{
	(my $dir) = shift; 
	opendir (DIR, $dir) || print "couldn't open $dir\n"; 
	my @dirAr = sort{lc($a) cmp lc($b)}(readdir(DIR));
	shift(@dirAr) while ($dirAr[0] =~ /^\./);
	close(DIR); 
	return @dirAr; 
}
