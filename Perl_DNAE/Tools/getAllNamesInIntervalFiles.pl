#Function: gets a unique list of all subfamilies in a directory full of interval files
#OUT: prints class ."\t". fam ."\t". subfam   to stdout. 
use strict; 

(my @invl_dirs) = @ARGV; 

my %subfams; 

#Get all names from all files
foreach my $dir(@invl_dirs){
	#get invl files from dir
	opendir (DIR, $dir) || die "$dir didn't open\n"; 
	my @invl_files = grep {/\.interval$/} readdir (DIR);
	closedir (DIR); 
	#get names from each file in dir
	foreach my $file(@invl_files){
		my $file_full = $dir ."/". $file; 
		open (FILE, $file_full) || die "open $file_full\n"; 
		while (my $line = <FILE>){
			next if $line =~ /^#/; 
			chomp $line; 
			my @fields = split(/\t/, $line); 
			$subfams{$fields[5]}{$fields[6]."\t".$fields[4]}=0; 
		}
		close(FILE); 
	}
}

#print unique subfamilies to stdout
foreach my $class (sort keys %subfams){
	foreach my $famSubfam(sort keys %{$subfams{$class}}){
		print $class ."\t". $famSubfam ."\n"; 
	}
}

