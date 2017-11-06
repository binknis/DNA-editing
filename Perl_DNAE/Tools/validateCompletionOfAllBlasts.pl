#checks that there are the same amount of blast queries in output file as sequences in the subfamilie's db file
use strict; 
use getAll; 
(my $org, my $class, my $fam, my $per_process) = @ARGV;
my $subfams = getAll::names($org, $class, $fam); 
die "no subfams for $org $class $fam\n" if $subfams == 0; 
my $home = $ENV{HOME}; 
my $out_dir = $home ."/". "validateBlastCompletion"; 
mkdir $out_dir; 
my $db_out = $out_dir ."/". "db_seq_count_per_fam_".$class.".txt"; 
my $blast_out = $out_dir ."/". "blast_output_count_per_fam_".$class.".txt"; 

if ($per_process){
	$db_out .= ".".$$; 
	$blast_out .= ".".$$; 
}

foreach my $subfam (@$subfams){
	###get amount of sequences in db file ###
	my $db_file = "$home/Data/$org/$class/db/files_".$fam."/Seq_".$subfam;
	$db_file =~ s/\(/\\\(/g; $db_file =~ s/\)/\\\)/g; 
	system ("echo -n $org $class $fam $subfam ' ' >> $db_out"); 
	system ("tac $db_file \| grep -m 1 '>' \| sed 's/>//' \| sed 's/=.*//' >> $db_out"); 
	
	### get amount of blasts completed for subfamily ### 
	my $zipped_blast = "$home/Data/$org/$class/results/blasts/$fam/Seq_".$subfam.".gz";
	$zipped_blast =~ s/\(/\\\(/g; $zipped_blast =~ s/\)/\\\)/g; 
	system ("echo -n $org $class $fam $subfam ' ' >> $blast_out"); 
	# system("gunzip -c $zipped_blast \| tac \| grep -m 1 Query= \| sed 's/Query= //' \| sed 's/=.*//' >> $blast_out"); 
	system("gunzip -c $zipped_blast \| grep -c Query= >> $blast_out"); 
}