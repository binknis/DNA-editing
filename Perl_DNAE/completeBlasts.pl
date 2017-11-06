#complete blasts for an organism which it's runClusterFinder.pl was stopped prematurely.
#Note: always does delIndexFiles!

#Note: MAKING THIS PARALLEL (= being able to run two processes on the same Class in parallel)

use strict; 
use File::Path qw(mkpath);
use File::Copy;
use IO::Compress::Gzip qw(gzip);
use IO::Uncompress::Gunzip qw(gunzip);
use lib "/home/alu/binknis/Perl_DNAE"; 
use removeIndexFiles; 

(my $organism, my $classFlag, my $delIndexFiles, my $cores) = @ARGV;

my @argClasses = @ARGV;
$cores = 1 if ($cores == 0); 

#Read Class names from file 
my $classDir = "../Data/" . $organism; 
opendir(CLASSES, $classDir) || print "classdir didn't open\n";
my @classList = sort{lc($a) cmp lc($b)}(readdir(CLASSES));
shift(@classList) while ($classList[0] =~ /^\./); #erase "." and ".." links
closedir(CLASSES);

#Retain only classes that are listed in the classFile (if classFlag == 1)
if ($classFlag){
        &getClasses(\@classList, $organism, \@argClasses); 
        if (scalar(@classList) == 0) { #No classes were retained - some kind of problem
                print "No classes retained in class-list for $organism\n"; 
        }
}

#erase all index files where clusters will be searched for (activate only when no other programs
#       are running BLAST on Names in this class). 
if ($delIndexFiles){
        removeIndexFiles::remove($organism, \@classList);
}

foreach my $class (@classList) {
        my $familyDir = $classDir ."/$class/db"; 
        opendir(FAMILY, $familyDir) || print "Family directory for $familyDir didn't open\n";
        my @familyList = sort{lc($a) cmp lc($b)}(readdir(FAMILY));
        shift(@familyList) while ($familyList[0] =~ /^\./); 
        closedir(FAMILY);


        foreach my $family (@familyList) {
                next unless (-d $familyDir ."/$family" ); #skip non-directory files
                my $nameDir = $familyDir . "/$family";
                next if ($family !~ /files_(\S+)/); #skip directories that don't contain sequences
                $family = $1; #erase "files_" prefix from family-name


                #create family's blast output files
                mkpath "../Data/$organism/$class/results/blasts/$family"; 


                opendir(NAME, $nameDir);
                my @nameList =  sort{lc($a) cmp lc($b)}(readdir(NAME));
                shift(@nameList) while ($nameList[0] =~ /^\./); 
                closedir(NAME);
                foreach my $name (@nameList) {
                        next if ($name =~ /\.n(hr|in|sd|si|sq|nd|ni|tm)$/); #skip blast index files (used becouse removeIndexFiles was commented)
                        my $blast_file = "../Data/$organism/$class/results/blasts/$family/$name";
                        my $blast_archive = $blast_file . ".gz"; 
                        next if (-e $blast_archive || -e $blast_file); #skip names which blasts' exist (the 2nd "-e" enables two processes of completeBlasts.pl to work on same Class simulataniously)
                        #create file now so that 2 processes don't start blast for the same Name. 
                        open(OUT, ">$blast_file") || print "open $blast_file failed\n"; 
                        close(OUT) || print "didn't close $blast_file\n"; 
                        #formatdb and BLAST
                        my $formatdb = "/private/apps/bin/formatdb -i $nameDir/$name -p F -o T"; 
						$formatdb =~ s/\(/\\\(/g; $formatdb =~ s/\)/\\\)/g; #add backslashes before "(" and ")" symbols. 
                        system($formatdb);
                        my $blast = "/private/apps/bin/blastall -p blastn -d $nameDir/$name -i $nameDir/$name -e 1e-50 -S 1 -F F -v 0 -a $cores > $blast_file";
						$blast =~ s/\(/\\\(/g; $blast =~ s/\)/\\\)/g; #add backslashes before "(" and ")" symbols. 
                        system($blast);
                        #zip blast output file
                        gzip $blast_file => $blast_archive or print "gzip failed for $blast_file \n";
                        unlink $blast_file or print "couldn't unlink $blast_file \n"; 

                }
        }
}

#for consistency - erase the index files that were created by this program
#removeIndexFiles::remove($organism, \@classList); #commented for parallel use


sub getClasses{
        my $classes_ref = shift; 

        my $organism = shift; 
        my $argClassList = shift; 
        my %classesToRetain = ();  

        if ($argClassList){ #specific classes were argumented to command line
                foreach my $class(@{$argClassList}){
                        $classesToRetain{$class} = 1; 
                }
        } 
        else{ #no argumented classes - use classes.txt class (this was the initial method, the "if" is the new feature using argumented classes in command line)
                open(my $classes_fh, "classes.txt") || print "Couldn't open classes.txt in runClusterFinder.pl - Running for all classes of $organism\n";

                while (<$classes_fh>){
                        foreach my $class (split){
                                $classesToRetain{$class} = 1; 

                        }
                }
                close($classes_fh); 

        }

        #create new classList (intersection of classes in classFile and in classesToRetain)
        my @newClassList = (); 
        foreach my $class (@{$classes_ref}){
                if (exists $classesToRetain{$class}){
                        push(@newClassList, $class); 

                }
        }
        @{$classes_ref} = @newClassList;
}
