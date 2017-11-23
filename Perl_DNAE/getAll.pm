#Package: Many utilities for scripts working with Org/Class/Fam directory tree
#Important note: the getDataDir() was retained for backward compatability, but it should not be used. 
#				 It will look for a dir named "Data" somewhere up the dir tree of where this file is present.
use strict; 
use Cwd; 
package getAll; 
use FindBin;



#pf = path flag: 1 = full path, 0 = no path. 

sub organisms{
	(my $dataDir, my $pf) = @_; 
	$dataDir = getDataDir() unless $dataDir; 
	my $orgs = myReadDir($dataDir, $pf); 
	return $orgs; 
}

#list of classes in an organism
sub classes{
	(my $dataDir, my $org) = @_; 
	$dataDir = getDataDir() unless $dataDir; 
	my $dir = $dataDir ."/". $org; 
	my $classes = myReadDir($dir); 
	return $classes;
}

#list of class files
sub class_files{
	(my $dataDir, my $org, my $pf) = @_; 
	$dataDir = getDataDir() unless $dataDir; 
	my $dir = $dataDir ."/". $org; 
	my $classes = myReadDir($dir, $pf); 
	return $classes;
}

#list of families 
sub families{
	(my $dataDir, my $org, my $class) = @_; 
	$dataDir = getDataDir() unless $dataDir; 
	my $dir = $dataDir ."/". $org ."/". $class ."/db";  
	my $family_dir_files = myReadDir($dir); 
	my @families; 
	eval{
		@families = grep {/files_[^\/]+$/} @$family_dir_files; 
	}; 
	if($@){
		return 0;
	}
	for(my $i=0;$i<=$#families;$i++){
			$families[$i] =~ s/^files_//;
	}
	return \@families; 
}

#list of families in "files_$family" format
sub family_files{
	(my $dataDir, my $org, my $class, my $pf) = @_; 
	$dataDir = getDataDir() unless $dataDir; 
	my $dir = $dataDir ."/". $org ."/". $class ."/db";  
	my $family_dir_files = myReadDir($dir, $pf); 
	my @families = grep {/files_[^\/]+$/} @$family_dir_files; 
	return \@families; 
}
#list of names
#Note: produces error if no names were retreived
sub names{
	(my $dataDir, my $org, my $class, my $family) = @_; 
	$dataDir = getDataDir() unless $dataDir; 
	my $dir = $dataDir ."/". $org ."/". $class ."/db/files_".$family;  
	my $names = myReadDir($dir);
	my @names_no_indexes = grep {!/\.n(hr|in|sd|si|sq|nd|ni|tm)$/} @$names; 
	for(my $i=0;$i<=$#names_no_indexes;$i++){
		$names_no_indexes[$i] =~ s/^Seq_//;
	}
	return \@names_no_indexes;
}
#list of name files
sub name_files{
	(my $dataDir, my $org, my $class, my $family, my $pf) = @_; 
	$dataDir = getDataDir() unless $dataDir; 
	my $dir = $dataDir ."/". $org ."/". $class ."/db/files_" . $family ;  
	my $name_files = myReadDir($dir, $pf); 
	my @names_no_index_files = grep {!/\.n(hr|in|sd|si|sq|nd|ni|tm)$/} @$name_files; 
	return \@names_no_index_files; 
}

#Get hash of all cluster files in tabular format (in subdirs for each pair of nucs; see @nuc_pairs)
#Each nuc-pair will either have a array of filenames returned (no path) or 0 if no files were found for the nuc-pair
sub tabClustersClass {
	(my $dataDir, my $org, my $class, my $pf, my $subdir) = @_; 
	$dataDir = getDataDir() unless $dataDir; 
	my $dir = $dataDir ."/". $org ."/". $class ."/results";  
	$dir .= "/$subdir" if $subdir; 
	my @nuc_pairs = ("GA", "CT", "GC", "GT", "CA", "TA");
	my %clustFiles = (); 
	my $dir_perPair; 
	my $file_count = 0; 
	foreach my $n (@nuc_pairs){
		$dir_perPair = $dir ."/". $n; 
		my $clusts_ref = myReadDir($dir_perPair, $pf); 
		# print @$clusts_ref, "\n"; #***
		$file_count += scalar(@$clusts_ref) unless $clusts_ref eq '0'; 
		$clustFiles{$n} = $clusts_ref; 
	}
	return 0 unless $file_count; #no files found
	return \%clustFiles; #files found
}
#Get hash of all cluster files in tabular format (in subdirs for each pair of nucs; see @nuc_pairs)
#Each nuc-pair will either have a array of filenames returned (no path) or 0 if no files were found for the nuc-pair
#*** Probably works but should be checked !!!
sub tabClustersFam {
	(my $dataDir, my $org, my $class, my $fam, my $pf, my $subdir) = @_; 
	$dataDir = getDataDir() unless $dataDir; 
	my $dir = $dataDir ."/". $org ."/". $class ."/results";  
	$dir .= "/$subdir" if $subdir; 
	my @nuc_pairs = ("GA", "CT", "GC", "GT", "CA", "TA");
	my %clustFiles = (); 
	my $dir_perPair; 
	my $file_count = 0; 
	foreach my $n (@nuc_pairs){
		$dir_perPair = $dir ."/". $n; 
		my $clusts_ref = myReadDir($dir_perPair, $pf); 
		if ($clusts_ref != 0){
			my $regex = "clusters_".$org."_".$class."_".$fam."_"; 
			my @clusts_per_fam = grep (/$regex/, @$clusts_ref); 
			$file_count += scalar(@clusts_per_fam);  
			$clustFiles{$n} = \@clusts_per_fam; 
		}
	}
	return 0 unless $file_count; #no files found
	return \%clustFiles; #files found
}

sub clustersClass{
	(my $dataDir, my $org, my $class, my $pf, my $subdir) = @_; 
	$dataDir = getDataDir() unless $dataDir; 
	my $dir = $dataDir ."/". $org ."/". $class ."/results";  
	$dir .= "/$subdir" if $subdir; 
	my $clusts_ref = myReadDir($dir, $pf); 
	return 0 unless $clusts_ref; #return 0 if class doesn't exist in organism
	$dir =~ s/\//\\\//g; #modify "dir". needed for regex if pf==1.
	my $regex = $pf ? "^".$dir."/clusters_" : "^clusters_"; 
	my @clusters = grep {/$regex/} @$clusts_ref; 
	return 0 unless $#clusters > -1; #return 0 if class doesn't exist in organism
	return \@clusters;
}

sub clustersClassAllOrgs{
	(my $dataDir, my $class, my $pf, my $subdir) = @_; 
	$dataDir = getDataDir() unless $dataDir; 
	my @clusts = (); 
	my $org_ref = organisms($dataDir); 
	foreach my $org (@$org_ref){
		my $clusts_ref = clustersClass($dataDir, $org, $class, $pf, $subdir);  
		next unless $clusts_ref; #organism doesn't have the class 
		push(@clusts, @$clusts_ref );
	}
	return \@clusts;
}

sub clustersFam{
	(my $dataDir, my $org, my $class, my $family, my $pf, my $subdir) = @_;
	$dataDir = getDataDir() unless $dataDir; 
	my $clusts_all = clustersClass($dataDir, $org, $class, $pf, $subdir); 
	return 0 unless $clusts_all; #return 0 if class isn't present in organism
	my @clusts_fam = (); 
	foreach my $cf (@$clusts_all){ #retain only cluster files of this family
		if ($cf =~ /clusters_([^_\/]+)_([^_\/]+)_([^\/]+)_(1e-\d+)_(\d+)(_control)?.txt$/){
			push (@clusts_fam, $cf) if $family eq $3; 
		}
	}
	return 0 unless $#clusts_fam > -1; #return 0 if family isn't present in class of organism
	return \@clusts_fam; 
}

sub clustersFamAllOrgs{
	(my $dataDir, my $class, my $family, my $pf, my $subdir) = @_;
	$dataDir = getDataDir() unless $dataDir; 
	my @clusts = (); 
	my $org_ref = organisms($dataDir); 
	foreach my $org (@$org_ref){
		my $clusts_ref = clustersFam($dataDir, $org, $class, $family, $pf, $subdir);  
		next unless $clusts_ref; #organism doesn't have the class 
		push(@clusts, @$clusts_ref);
	}
	return \@clusts;
}

# sub clustersParamsOrg{
	# (my $org, $class, )
# }

# sub clustersParamsAll{
	# (my $org, $class, )
# }


sub clustStats{
	(my $dataDir, my $org, my $class, my $pf, my $subdir) = @_; 
	$dataDir = getDataDir() unless $dataDir; 
	my $dir = $dataDir ."/". $org ."/". $class ."/results";  
	$dir .= "/$subdir" if $subdir; 
	my $results_ref = myReadDir($dir, $pf); 
	return 0 unless $results_ref; #return 0 if class doesn't exist in organism
	$dir =~ s/\//\\\//g; #modify "dir". needed for regex if pf==1.
	my $regex = $pf ? "^".$dir."/cluster_stats_" : "^cluster_stats_"; 
	my @clusterStats = grep {/$regex/} @$results_ref; 
	return \@clusterStats;
}

# sub clustStatsAllOrgs{
	
# }

#cluster-stats files of a Fam in a specific Organim's Class
sub clustStatsFam{
	(my $dataDir, my $org, my $class, my $family, my $pf, my $subdir) = @_;
	$dataDir = getDataDir() unless $dataDir; 
	my $clustStats = clustStats($dataDir, $org, $class, $pf, $subdir); 
	return 0 unless $clustStats; #return 0 if no clustStats files were found
	# my $regex = $class."_".$family.".txt\$"; 
	my $regex = $class."_".$family.".txt\$"; 
	my @clusterStatsFam = grep {/$regex/} @$clustStats; 
	return 0 unless $#clusterStatsFam > -1; 
	return \@clusterStatsFam;
}

sub clustStatsFamAsFilename{
	(my $dataDir, my $org, my $class, my $family, my $pf, my $subdir) = @_;
	$dataDir = getDataDir() unless $dataDir; 
	my $clustStats = clustStats($dataDir, $org, $class, $pf, $subdir); 
	return 0 unless $clustStats; #return 0 if no clustStats files were found
	# my $regex = $class."_".$family.".txt\$"; 
	my $regex = $class."_".$family.".txt\$"; 
	my @clusterStatsFam = grep {/$regex/} @$clustStats; 
	return 0 unless $#clusterStatsFam > -1; 
	return $clusterStatsFam[0];
}

sub clustStatsFamAllOrgs{
	(my $dataDir, my $class, my $family, my $pf, my $subdir) = @_;
	$dataDir = getDataDir() unless $dataDir; 
	my @clustStats = ();
	my $org_ref = organisms($dataDir); 
	foreach my $org (@$org_ref){
		my $clustStats_ref = clustStatsFam($dataDir, $org, $class, $family, $pf, $subdir);  
		next unless $clustStats_ref; #organism doesn't have the class 
		push(@clustStats, @$clustStats_ref );
	}
	return \@clustStats;
}

#Returns a hash-ref with nucleotide frequencies by subfamily in the specified family
sub nucFreqFamBySubfam{
	(my $dataDir, my $org, my $class, my $family) = @_; 
	$dataDir = getDataDir() unless $dataDir; 
	my $nuc_fam_file = $dataDir ."/".$org."/".$class."/db/Nuc_".$family.".txt"; 
	my $nuc_lines_ref = lines($nuc_fam_file); 
	my %nuc_by_subfam = (); 
	foreach my $line (@$nuc_lines_ref){
		(my $subfam) = ($line =~ /(\S+) /);
		my $temp = $';
		my @arr = split(/ /,$temp);
		@arr = @arr[0..3];
		$nuc_by_subfam{$subfam} = [@arr];
	}
	return \%nuc_by_subfam; 
}

#input: $org, $class, $family, $subfam
#returns: an array-ref with the nucleotide frequencies of the subfam
sub nucFreqSubfam{
	(my $subfam) = pop(@_);  
	my $nuc_freq_fam = nucFreqFamBySubfam(@_);  
	return $nuc_freq_fam->{$subfam};
}
#Returns: an array-ref with the nucleotide frequencies of the family (needs the external file created by calcFamNucFreq.pl)
sub nucFreqFam{
	(my $org, my $class, my $family) = @_; 
	my $nuc_by_fam_file = "/home/alu/binknis/DNA_editing_results/db_stats/Family_nuc_freq_all_orgs.txt"; 
	print "Warning: running nucFreqFam with $nuc_by_fam_file - make sure it's updated!"; #***
	my $fam_lines = lines($nuc_by_fam_file); 
	my $regex = $org ."\t". $class ."\t". $family;
	(my $fam_nuc_line) = grep (/$regex/ , @$fam_lines); 
	my $regex2 = $org."\t".$class."\t".$family."\s+"; 
	$fam_nuc_line =~ /$regex2/;
	(my @freqs) = ( $' =~ /(\d+)/g );
	return \@freqs;
}

#Read all lines of a file into an array and return its ref.
sub lines{
	(my $file_full_path, my $chomp) = @_; 
	open (FILE, "<$file_full_path") || return 0; 
	my @lines = <FILE>;
	close(FILE); 
	if ($chomp){
		for (my $i=0; $i<$#lines; $i++){
			chomp $lines[$i];
		}
	}
	return \@lines; 
}

#get the "Data" directory. created so that you don't have to be one level under Data to call this script's functions
sub getDataDir{
	my $perlDir = "$FindBin::Bin"; 
	my $dataDir = $perlDir; 
	while(not -d "$dataDir/Data"){
		($dataDir) = ($dataDir =~ /(\S+)\/\S+/); 
	}
	return "$dataDir/Data";
}

#Function: returns a list of files from a directory (excluding "." and ".." links).
#input: full path to directory and path_flag: 1 - add path to file, 0 - file w/o path
#ret: ref to array of file names (full path, if path_flag = 1)
#		0 - if dir doesn't exist
sub myReadDir{
	(my $dir, my $path_flag)  = @_; 
	return 0 unless -d $dir; #dir doesn't exist
	opendir(DIR, $dir) || print "$dir didn't open\n";
	my @files = sort{lc($a) cmp lc($b)} grep {!/^\./} readdir(DIR); #erase "." and ".." links
	closedir(DIR);
	if ($path_flag){
		for (my $i=0;$i<=$#files;$i++){ #add paths to filenames
			$files[$i] = $dir ."/". $files[$i]; 
		}
	}
	return \@files; 
}

		
1; 