#Function: erases all files (db and results) of a family or names in one family
#			For family - deletes all files completely. 
#			For names - deletes only relevent lines in the files (db and results)
#Note: doesn't erase names from files or files created by advanced analysis (only original db files, cluster files and cluster_stat files).
#for Name deletion - you MUST create the cluster_stat files again.

use strict; 
use getAll; 
use File::Copy; 
(my $org, my $class, my $fam) = splice(@ARGV, 0, 3); 
(my @names) = @ARGV; 

my $dataDir = getAll::getDataDir(); 
my $dbDir = $dataDir ."/". $org ."/". $class ."/db"; 
my $lenFile = $dbDir ."/Len_".$fam.".txt"; 
my $lenStatsFile = $dbDir ."/LenStats_".$fam.".txt"; 
my $nucFile = $dbDir ."/Nuc_".$fam.".txt"; 

my $resultDir = $dataDir ."/". $org ."/". $class ."/results"; 

unless (@names) {#### no names specified - delete whole FAMILY!! ####
	#erase fam db files
	my $names_ref = getAll::names($org, $class, $fam); #save names for checkCount deletion (if exist)
	system ("rm -r $dbDir/files_".$fam); 
	foreach my $dbFile ($lenFile, $lenStatsFile, $nucFile){
		system("rm $dbFile");
	}
	#erase fam cluster files
	system("rm ".$resultDir."/cluster*_".$org."_".$class."_".$fam."[_\.]*"); 
	
	#erase blast files
	system("rm -r $resultDir/blasts/$fam"); 
	
	#erase checkCount files (if they exist)
	foreach my $name (@$names_ref){
		system("rm -f $resultDir/checkCount/Seq_$name"); 
	}
}
else{ #### NAMES specified - del specific names in family ####
	my %namesHash = (); 
	foreach my $name (@names){
		$namesHash{$name}=1;	
		#erase sequence files
		system ("rm $dbDir/files_".$fam."/Seq_".$name);
		#erase names' blast files
		system("rm $resultDir/blasts/$fam/Seq_".$name.".gz");
		
		#erase checkCount files (if they exist)
		foreach my $name (@names){
			system("rm -f $resultDir/checkCount/Seq_$name"); 
		}
	}
	
	#erase names db files
	foreach my $dbFile ($lenFile, $lenStatsFile, $nucFile){
		eraseNamesLinesFromFile($dbFile, \%namesHash); 
	}		
	#erase names' clusters from cluster files
	my $clusterFiles_ref = getAll::clustersFam($org, $class, $fam, 1); 
	if ($clusterFiles_ref ne '0'){
		foreach my $cf (@$clusterFiles_ref){
			eraseNamesClustersFromCF($cf, \%namesHash); 
		}
	}
	else {
		print "no cluster files retreived\n"; 
	}
}

########################
##### SUBROUTINES ######
########################

#Input: 1. a filename 2. hash of names to erase their lines from the file
#erases lines with the name in the head of the line (and then a whitespace)
sub eraseNamesLinesFromFile{
	(my $fileName, my $namesHash_ref) = @_; 
	my $orig = $fileName.".orig"; 
	move($fileName, $orig); 
	open(FINAL, ">$fileName") or die "couldn't open $fileName\n"; 
	open(ORIGINAL, "<$orig") or die "couldn't open $orig\n"; 
	while (my $line = <ORIGINAL>){
		next unless $line =~ /^(\S+)/; 
		my $name = $1; 
		unless (exists $namesHash_ref->{$name}){
			chomp $line; 
			print FINAL $line ."\n"; 
		}
	}
	close(ORIGINAL); 
	close(FINAL); 
	unlink($orig); 
}


#Function: erases all clusters generated from names in the name hash from a clusters file
#Input: 1. cluster-file filename. 2. hash of names to delete their clusters. 
sub eraseNamesClustersFromCF{
	(my $cfName, my $namesHash_ref) = @_; 
	my $orig = $cfName.".orig"; 
	move($cfName, $orig); 
	open(FINAL, ">$cfName") or die "couldn't open $cfName\n"; 
	open(ORIGINAL, "<$orig") or die "couldn't open $orig\n"; 
	while (my $line = <ORIGINAL>){
		if ($line =~ /^Found/){ #for each cluster
			my $cluster = $line; 
			do {
				$line = <ORIGINAL>; 
				$cluster .= $line;
			}while ($line !~ /^End cluster/);
			#print the cluster only if its name isn't in the name's list
			(my $name) = $cluster =~ /name = \S+=(\S+)/; 
			print FINAL $cluster unless exists $namesHash_ref->{$name}; 
		}
	}
	close(ORIGINAL); 
	close(FINAL); 
	unlink($orig); 
}






