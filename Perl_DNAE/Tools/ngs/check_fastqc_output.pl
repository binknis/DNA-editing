#!/private/apps/bin/perl516 -w
use strict;
use warnings;
#Check if the quality test of fastqc passed

die "$0 fastqc_text_file \n" unless @ARGV;
my $input_file = $ARGV[0];
open(IN,$input_file);

my $line;
my @l;

while ($line = <IN>)
{
	chomp($line);
	@l = split("\t",$line);
	if($l[1] eq "Per base sequence quality" && $l[0] eq "FAIL")
	{
		warn ("$input_file Per base sequence quality test failed");
		die;
	}
	if($l[1] eq "Per sequence quality scores" && $l[0] eq "FAIL")
	{ 
		warn ("$input_file Per sequence quality scores test failed");
		die;
	}
}
close(IN);

