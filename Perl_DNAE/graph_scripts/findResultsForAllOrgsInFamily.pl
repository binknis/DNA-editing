use strict; 

(my $class, my $family, my $nonDefaultPath, my @orgs) = @ARGV; 

my $dir = "../Data"; 
$dir = $nonDefaultPath unless $nonDefaultPath eq ""; 

my %bestParamData = (); 
my %paramsToOrgs = ();

#open clusterStats and find existing parameters
foreach my $org (@orgs){
	my $resultDir = "$dir/$org/$class/results"; 
	opendir(RESULTS, $resultDir) || die "couldn't open $resultDir\n"; 
	readdir(); #***
	close(RESULTS); 
	
	

}

opendir(CLASSES, $classDir) || print "classdir didn't open\n";
@classList = sort{lc($a) cmp lc($b)}(readdir(CLASSES));
shift(@classList) while ($classList[0] =~ /^\./); #erase "." and ".." links
closedir(CLASSES);