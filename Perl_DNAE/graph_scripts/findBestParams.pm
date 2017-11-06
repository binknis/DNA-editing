use strict;
(my $organism , my $class ) = @ARGV;
my $dir = "../Data/". $organism . "/" .$class ."/results";
my $line;
my @arr;
my $p_value;
my $th;
my $now=0;
my $next=0;
my %hash=();
open (FILE_STAT ,"$dir/cluster_stats.txt")  || die("the cluster stats is error");
open (my $best_file,">>$dir/bestParamsByFamily.txt") || die("the best file is error");
while($line = <FILE_STAT> )
{
 @arr = split(/\t/,$line);
 if ($arr[4]== 0 && $arr[5] == 0)
    {next;}
 if (not exists $hash{$arr[0]})
   {
	$hash{$arr[0]} = $arr[1]." ".$arr[2]." ".$arr[3];
    $now = abs (10 - $arr[3]);
   }
 else
  {
    $next = abs (10 - $arr[3]);
    if($next < $now)
	{
	 $hash{$arr[0]} = $arr[1]." ".$arr[2]." ".$arr[3];
     $now = abs (10 - $arr[3]);	 
	}	
  }
  
  

}
 while ( my ($key, $value) = each(%hash) ) {
        print $best_file "$key: $value\n";
    }

 
close (FILE_STAT);
close($best_file);