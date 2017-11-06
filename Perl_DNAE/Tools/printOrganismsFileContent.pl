 #ARG: name of organism (reached by: ../Data/$Organism)
 #function: prints file tree of the organism ($class/$family/$name)
 #INPUT: name of organism which has been sorted by sortGenome.pl thus creating its file-tree
 
 use strict; 
 
(my $organism, my $pvalue, my $th) = @ARGV;
 
###LOOP 1: each class in organism###
my $classDir = "../Data/" . $organism; 
opendir(CLASSES, $classDir) || print "classdir didn't open\n";
my @classList = sort{lc($a) cmp lc($b)}(readdir(CLASSES));
closedir(CLASSES);
foreach my $class (@classList) {
	next if ($class =~ /^\./); #skip "." and ".." links

	print "$class: \n"; 
	###LOOP 2: each family in class###
	my $familyDir = $classDir ."/$class/db"; 
	opendir(FAMILY, $familyDir) || print "familydir didn't open\n";
	my @familyList = sort{lc($a) cmp lc($b)}(readdir(FAMILY));
	closedir(FAMILY);
	foreach my $family (@familyList) {
		next if ($family =~ /^\./); #skip "." and ".." links
		next unless (-d $familyDir ."/$family" ); #skip non-directory files

		my $nameDir = $familyDir . "/$family";
		
		next if ($family !~ /files_(\S+)/); #skip directories that don't contain sequences
		print "    ->$family\n"; 
	
		
		###LOOP 3: each name in family###
		opendir(NAME, $nameDir);
		my @nameList =  sort{lc($a) cmp lc($b)}(readdir(NAME));
		closedir(NAME);
		foreach my $name (@nameList) {
			next if ($name =~ /^\./); #skip "." and ".." links
			print "     ->      ->$name\n"; 
		
		}
		
	}
}

