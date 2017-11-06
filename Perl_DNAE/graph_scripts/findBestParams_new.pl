use strict;
(my $organism , my $class ) = @ARGV;
my $dir = "../Data/". $organism . "/" .$class ."/results";
my $line;
my @arr;
my $p_value;
my $th;
my $now=0;
my $next=0;
my @max=();
my %hash_res=();
my $value=0;
my $value1=0;
my $value2=0;
my $value3=0;
my %hash=();
my @results=();
open (FILE_STAT ,"$dir/cluster_stats.txt")  || die("the cluster stats is error");
open (my $best_file,">>$dir/bestParamsByFamily.txt") || die("the best file is error");
while($line = <FILE_STAT> )
{
 
 @arr = split(/\t/,$line);
 if ($arr[4]== 0 && $arr[5] == 0)
    {next;}
 
 
 if (not exists $hash{$arr[0]})
   {
	@results = ();
	for (my $i=0 ; $i<4; $i++)
	  {
	   $results[$i] = "Not found";
	  }
	$hash{$arr[0]} = $arr[1]." ".$arr[2]." ".$arr[3];
    $now = abs (10 - $arr[3]);
	$value=0;
    $value1=0;
    $value2=0;
    $value3=0;
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
   #Max(c4-c5) > 0
 $max[0] = $arr[4] - $arr[5];
 if($max[0] > $value)
   {
   $value=$max[0];
   $results[0] = $arr[1]." ".$arr[2]." ".$arr[3]; 
   } 
 #Max(c4) AND c5==0
 $max[1] = $arr[4];
 if($max[1] > $value1 && $arr[5]== 0)
   {
   $value1=$max[1];
   $results[1] = $arr[1]." ".$arr[2]." ".$arr[3]; 
   }
 #Max(c5-c4) > 0
 $max[2] = $arr[5] - $arr[4];
  if($max[2] > $value2)
   {
   $value2=$max[2];
   $results[2] = $arr[1]." ".$arr[2]." ".$arr[3]; 
   }
 #Max(c5) AND c4==0
 $max[3] = $arr[5];
 if($max[3] > $value3 && $arr[4]== 0)
   {
   $value3=$max[3];
   $results[3] = $arr[1]." ".$arr[2]." ".$arr[3]; 
   }
   
   
  $hash_res{$arr[0]} = \@results;

}
 while ( my ($key, $value) = each(%hash) ) {
        print $best_file "Family $key:\n";
		print $best_file "best FP: $value\n";
	    foreach my $parameters (@{$hash_res{$key}}){
			print $best_file "$parameters\n";
		}
    }

 
close (FILE_STAT);
close($best_file);