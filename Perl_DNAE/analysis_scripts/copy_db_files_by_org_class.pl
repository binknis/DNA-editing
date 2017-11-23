#Function: gets a organism and class name and copies its db files to a united file for all organisms and classes
#script to use to unite from all organisms:  
#	perl foreachOrgClassExec.pl "perl analysis_scripts/copy_db_files_by_org_class.pl xxx yyy"
use strict; 
use getAll; 

(my $dataDir, my $org, my $class) = @ARGV; 
my $fams = getAll::families($dataDir, $org, $class); 
die "No fams for: $org $class\n" if $fams eq 0; 

### copy file content ###
my $dbdir = "$dataDir/Data/$org/$class/db";
my $outdir = "/home/alu/binknis/DNA_editing_results/db_stats"; 
foreach my $fam (@$fams){
	#length file
	my $len = $dbdir ."/Len_".$fam.".txt"; 
	my $lenOut = $outdir ."/Len_all_orgs_classes.txt"; 
	printWithOrgClass($len, $lenOut, $org, $class, $fam); 
	#length stats file
	my $lenStats = $dbdir ."/LenStats_".$fam.".txt"; 
	my $lenStatsOut = $outdir ."/LenStats_all_orgs_classes.txt"; 
	printWithOrgClass($lenStats, $lenStatsOut, $org, $class, $fam); 
	my $nuc = $dbdir ."/Nuc_".$fam.".txt";
	my $nucOut = $outdir ."/Nuc_all_orgs_classes.txt"; 
	printWithOrgClass($nuc, $nucOut, $org, $class, $fam); 
}

sub printWithOrgClass { 
	(my $in, my $out, my $org, my $class, my $fam) = @_;
	open (IN, "<$in") or die "$in\n"; 
	open (OUT, ">>$out") or die "$out\n"; 
	
	my @lines = <IN>;
	foreach my $line (@lines){
		next if $line =~ /^#/; #***make sure it's correct
		$line =~ s/(\S+)\s//; #erase subfamily name from file 
		print OUT $org ."\t". $class ."\t". $fam ."\t". $1 . "\t" . $line; 
	}
	close(IN); 
	close(OUT); 
}

