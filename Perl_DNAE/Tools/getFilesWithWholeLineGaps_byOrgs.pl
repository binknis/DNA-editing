#Find files that are problematic for SearchIO module .
#problematic line: "Query    --------------------------------  " (Query with gaps only. doens't have numbers at head and end of file)
#NOTE: use this by REDIRECTING the output to file
#NOTE: subsequent manual removal of problematic HSPs is required!
#input: 1. The directory of the organisms; 2. a list of organisms to search in for problematic BLAST files.  
use strict; 

(my $orgDir, my @orgInput) = @ARGV; 
my $orgRegex = "(" . join ("|",@orgInput) . ")"; 


opendir(ORGS, $orgDir) || print "$orgDir didn't open\n"; 
my @organisms = sort{lc($a) cmp lc($b)}(readdir(ORGS));
shift(@organisms) while ($organisms[0] =~ /^\./); 
closedir(ORGS);

#organism
foreach my $org (@organisms){
	next unless $org =~ /$orgRegex/;
	my $classDir = $orgDir ."/". $org;
	my @classes = &getFilesInDir($classDir); 
	#class
	foreach my $class (@classes){
		next unless $class =~ /SINE/;
		my $familyDir = $classDir ."/". $class . "/results/blasts"; 
		my @fams = &getFilesInDir($familyDir);
		# family
		foreach my $files_fam (@fams){
			my $nameDir = $familyDir."/".$files_fam;
			my @seq_names = &getFilesInDir($nameDir);
			foreach my $seq_name (@seq_names){
				my $full_seq_name = $nameDir . "/" . $seq_name; 
				$full_seq_name =~ s/\(/\\\(/g; $full_seq_name =~ s/\)/\\\)/g; #add backslashes before "(" and ")" symbols. 
				print $full_seq_name . " "; 
				my $grep = "gunzip -c " . $full_seq_name . ' | grep -c -P "^Query\s+[-]+\s*$" ';
				system ($grep); 
			}
		}
	}
}

sub getFilesInDir{
	(my $dir) = shift; 
	opendir (DIR, $dir) || print "couldn't open $dir\n"; 
	my @dirAr = sort{lc($a) cmp lc($b)}(readdir(DIR));
	shift(@dirAr) while ($dirAr[0] =~ /^\./);
	close(DIR); 
	return @dirAr; 
}