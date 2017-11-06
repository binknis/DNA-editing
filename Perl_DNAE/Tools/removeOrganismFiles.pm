####   Erase subroutine   ####
#Function: Erases all files made by formatdb from a specified genome.
#"classes" arg: if passed - deletes only from those classes. 
# 				if not - deletes from all classes in Organism's directory. 
#genome file-tree format: sorted by sortGenome.pl 
#suffixes are listed in %sfxToDel. 
package removeIndexFiles;
use strict; 

sub removeIndexFiles{
 
	my $organism = shift;
	my $classes = shift; 
	my @classList; 
	#create hash of filenames suffixes to delete:
	my %sfxToDel = ('.nhr' => 1, '.nin' => 1, '.nsd'=> 1, '.nsi'=> 1, '.nsq'=> 1); 
	 
	###LOOP 1: each class in organism###
	my $classDir = "../Data/" . $organism; 

	if ($classes){ #use class list passed as argument. 
		@classList = @{$classes}; 
	}
	else{ #no classList argument => extract class data from Organism's directory.
		opendir(CLASSES, $classDir) || print "classdir didn't open\n";
		@classList = sort{lc($a) cmp lc($b)}(readdir(CLASSES));
		shift(@classList) while ($classList[0] =~ /^\./); #erase "." and ".." links
		closedir(CLASSES);
	}

	foreach my $class (@classList) {
		###LOOP 2: each family in class###
		my $familyDir = $classDir ."/$class/db"; 
		opendir(FAMILY, $familyDir);
		my @familyList = sort{lc($a) cmp lc($b)}(readdir(FAMILY));
		shift(@familyList) while ($familyList[0] =~ /^\./); 
		closedir(FAMILY);
		foreach my $family (@familyList) {
			next unless (-d $familyDir ."/$family" ); #skip non-directory files
			next if ($family !~ /files_(\S+)/); #skip directories that don't contain sequences
			my $nameDir = $familyDir . "/$family";
			###LOOP 3: each name in family###
			opendir(NAME, $nameDir);
			my @nameList =  sort{lc($a) cmp lc($b)}(readdir(NAME));
			shift(@nameList) while ($nameList[0] =~ /^\./); 
			closedir(NAME);
			foreach my $name (@nameList) {
				my $end = length($name)-1; 
				if (exists $sfxToDel{substr($name,$end-3,$end)}){
					unlink($nameDir . "/$name"); 
				}
			}
			
		}
	}

}

sub removeResults{
	my $organism = shift; 
	my $classes = shift; 
	my $flag = ;  #***
	
	###LOOP 1: each class in organism###
	my $classDir = "../Data/" . $organism; 

	if ($classes){ #use class list passed as argument. 
		@classList = @{$classes}; 
	}
	else{ #no classList argument => extract class data from Organism's directory.
		opendir(CLASSES, $classDir) || print "classdir didn't open\n";
		@classList = sort{lc($a) cmp lc($b)}(readdir(CLASSES));
		shift(@classList) while ($classList[0] =~ /^\./); #erase "." and ".." links
		closedir(CLASSES);
	}
	
	foreach my $class (@classList) {
		###LOOP 2: each family in class###
		my $resultsDir = $classDir ."/$class/results"; 
		opendir(RESULTS, $resultsDir);
		my @resultsList = readdir(RESULTS);
		shift(@resultsList) while ($resultsList[0] =~ /^\./); 
		closedir(RESULTS);
		foreach my $result (@resultsList) {
				
			#*** decide what to delete 
			
			unlink(); 
		}
	}
	
	
}


1;