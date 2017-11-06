####   Erase subroutine   ####
#Function: Erases all files made by AnalyzeBlastByLength.pl (cluster and checkCount files)
#"classes" arg: if passed - deletes only from those classes. 
# 				if not - deletes from all classes in Organism's directory. 
#erases ALL non-directory files in results!!!
use strict; 
use File::Path;

my $organism = shift;
my @classList = @ARGV;
 
###LOOP 1: each class in organism###
my $classDir = "../Data/" . $organism; 

if ($#classList == -1){ #no classList argument => extract class data from Organism's directory.
	opendir(CLASSES, $classDir) || print "$classDir didn't open\n";
	@classList = sort{lc($a) cmp lc($b)}(readdir(CLASSES));
	shift(@classList) while ($classList[0] =~ /^\./); #erase "." and ".." links
	closedir(CLASSES);
}

foreach my $class (@classList) {
	my $clustDir = $classDir ."/$class/results"; 
	opendir(CLUST, $clustDir) || print "$clustDir didn't open!\n";
	my @clustList =  sort{lc($a) cmp lc($b)}(readdir(CLUST));
	shift(@clustList) while ($clustList[0] =~ /^\./); 
	closedir(CLUST);
	foreach my $file (@clustList) { #each cluster file
		#skip directories that don't contain sequences
		my $fileFullPath = "$clustDir/$file"; 
		if (-f $fileFullPath){ #erase 
			unlink($fileFullPath); 
		}
		elsif ($file =~ /^checkCount/ && -d $fileFullPath){
			rmtree($fileFullPath); 
		}
	}
		
}
