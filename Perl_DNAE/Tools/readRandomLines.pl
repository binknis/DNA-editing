#Function: get random lines from a file
#Input: $file - input file
#		$num - number of lines to get
#		$percent - flag if $num is percent (1) or absolute number (0)
#Notes: 1. Random numbers are chosent without replacement
#		2. if $num exceeds num_of_lines_in_file then $num will be reset to num_of_lines_in_file (i)
#		3. For percentage - when calculating number of lines the digits after decimal are truncated; Min 1 lines will always be printed. 
#		4. Doesn't read from stdinput
use strict; 
use List::Util qw(max);

(my $file, my $num, my $percent) = @ARGV; 
#check input
if ($num==0){
	die "bad input - num must be greater than 0\n"; 
}
if ($percent and ($num>100)){
	die "bad input - percent can't exceed 100\n"; 
}

#Get number of lines in file
open (FILE, $file) or die "open $file\n"; 
while (<FILE>){}
my $numLines = $.; 
close(FILE); #must be after use of $. (otherwise the value is reset)
if ($percent){
	$num = max(1,int($numLines * $num / 100)); #max is to avoid zero
}
if ($num > $numLines){ #num of random lines can't exceed the number of lines in file
	$num = $numLines; 
}
  
#generate unique random numbers
my $num_generated = 0; 
my %nums = (); 
my $randNum; 
while ($num_generated < $num){
	$randNum = int(rand($numLines)) + 1; #+1 for min 1 max
	unless(exists $nums{$randNum}){ #no replacement
		$nums{$randNum}=1; 
		$num_generated++; 
	}
}

#print the random lines to stdout
open (FILE, $file) or die "open $file\n"; 
while (my $l = <FILE>){
	if (exists $nums{$.}){
		print $l; 
	}
}
