#INPUT: Organism, class, family, name. 
#OUTPUT: amount of queries done out of the whole Name file. 
use strict; 


(my $organism, my $class, my $family, my $name) = @ARGV; 

my $seqFile = "../Data/$organism/$class/db/files_$family/Seq_$name";
my $blastFile = "../Data/$organism/$class/results/blasts/$family/Seq_$name";

print "Total sequences for $name: "; 
system("grep -c \">\" $seqFile"); 

print "Total blast complete: "; 
system("grep -c ^Query= $blastFile"); 



