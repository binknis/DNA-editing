#Function: Find duplicate coordinates in same file (Typically with different annotation, which is what causes the redundance)
#whole_dir_flag: if to run for all interval files in the directory specified
#Intreval files must have format "\S+.interval"
####FINISHED - NEED TO CHECK THAT IT WORKS
use strict; 
(my $interval_file_or_dir, my $whole_dir_flag) = @ARGV; 
my %coords = (); 
my @files; 
if ($whole_dir_flag){ #run for all interval files in specified dir
	opendir (INVL_DIR, $interval_file_or_dir) | die "open $interval_file_or_dir\n";  
	@files = grep {/\.interval$/} readdir (INVL_DIR); 
	close(INVL_DIR); 
}
else{ #run for specified interval file
	push (@files, $interval_file_or_dir); 
}

foreach my $interval_file (@files){
	(my $out_file) = $interval_file =~ /^(\S+)\.interval$/; 
	open (INVL, $interval_file) || die "open $interval_file\n"; 
	 while(my $line = <INVL>){
		chomp $line; 
		my @fields = split (/\s+/,$line); 
		my $coord = join ('_', @fields[0 .. 3]); 
		$coords{$coord}=0 unless exists $coords{$coord}; 
		$coords{$coord}++; 
	}
	close(INVL);

	open (INVL, $interval_file) || die "open $interval_file\n"; 
	open (OUT, ">$out_file") || die "open $out_file\n"; 
	 while(my $line = <INVL>){
		chomp $line; 
		my @fields = split (/\t/,$line); 
		my $coord = join ('_', @fields[0 .. 3]); 
		if ($coords{$coord} > 1){
			print OUT $line ."\n"; 
		}
	}
	close(OUT); 
	close(INVL);
}