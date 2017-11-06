#Function: deletes redundant lines from db files (LenStats, Len, Nuc). 
#			Saves the last listing of each subfamily. 
use strict; 
use getAll;
use File::Copy;  
(my $org, my $class) = @ARGV; 
my @prefixes = ("LenStats", "Len", "Nuc");
my $fams = getAll::families($org, $class); 
for my $pref(@prefixes){
	for my $fam(@$fams){
		my %line_by_subfam = (); 
		my $file = "../Data/$org/$class/db/".$pref ."_". $fam .".txt";
		my $lines = getAll::lines($file);
		for my $line(@$lines){
			chomp $line; 
			(my $subfam) = $line =~ /^\s*(\S+)/;
			$line_by_subfam{$subfam} = $line if ($subfam =~ /\S+/); 
		}
		my $temp_file = $file .".temp";
		open (TEMP, ">$temp_file") || die "couldn't open temp file: $file.temp\n"; 
		for my $sf (sort keys %line_by_subfam){
			print TEMP $line_by_subfam{$sf} ."\n"; 
		}
		close TEMP;
		move ($temp_file, $file);
	}
}