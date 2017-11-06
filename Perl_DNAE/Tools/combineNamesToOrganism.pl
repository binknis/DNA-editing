#For a family created seperately after running RepeatMasker
#Function: moves all data found in the family to the actual organism's directory
#notes: 1. moves only raw clusters (no subdirs)
#		2. If same Family and Subfamily exist in new_org
#			a. Seq files will be replaced. 
#			b. Redundancy will be found in the 3 db files (2 different lines for 1 subfamily) -- NEED TO DELETE MANUALLY
#			c. Blasts will be replaced
#			d. Checkcount files will be replaced
#		3. If same Family exists - CLUSTER files will be appended. Otherwise, created.
#		4. Notifies (stdout) if cluster files don't exist for the fam
#		5. PROBLEM if name needs to be changed and cluster files exist- doesn't change annotation IN cluster files! names must have 'family' changed per cluster
use strict; 
use File::Path qw(mkpath);
use getAll; 
(my $org, my $class, my $fam, my $new_org, my $new_fam, my @names) = @ARGV; 

my %namesH = map { $_ => 1 } @names;
my $old_dir = "../Data/$org/$class"; 
my $new_dir = "../Data/$new_org/$class";

### copy db files ###
mkpath $new_dir."/db"; 
#1. seq files
#get all Seq files from old fam lib
my $seqDir  = $old_dir ."/db/files_".$fam; 
opendir(SEQDIR, $seqDir) or die "didn't open $seqDir for $fam\n"; 
my @seqFilesAll = grep {!/\./} readdir (SEQDIR);  #only names without any dot characters. discards index files
closedir(SEQDIR);

#retain only seq files of names wanted
my @seqFiles = (); 
for my $sf(@seqFilesAll){
 $sf =~ /Seq_(\S+)/; 	
 push (@seqFiles, $sf) if exists $namesH{$1}; 
}
die "no names retained for $org $class $fam $new_org $new_fam - the names specified probably don't exist in the family specified\n" if $#seqFiles == -1; 

my $new_seqDir = $new_dir ."/db/files_".$new_fam;
mkpath $new_seqDir; 
foreach my $seqFile (@seqFiles){
	open (SEQ, "<".$seqDir ."/". $seqFile) or die "didn't open seqDir for $org $fam\n"; 
	open (NEWSEQ, ">".$new_seqDir ."/". $seqFile) or die "didn't open new seqDir for $org $fam\n"; 
	while (my $line = <SEQ>){
		chomp $line; 
		if ($line =~ /^>/){
			(my $head, my $tail) = $line =~ /^(>\S+=)[^=]+(=[^=]+)$/; 
			$line = $head .$new_fam. $tail;
		}
		print NEWSEQ $line ."\n";
	}
	close (SEQ); 
	close(NEWSEQ); 
}

#2. Len, LenStats, Nuc
my @prefixes = ("LenStats", "Len", "Nuc"); 
foreach my $pref (@prefixes){
	my $file = $old_dir."/db/$pref"."_".$fam.".txt"; 
	my $new_file = $new_dir."/db/$pref"."_".$new_fam.".txt"; 
	my $lines = getAll::lines($file); 
	for my $line (@$lines){
		chomp $line; 
		$line =~ /^(\S+)/;
		system("echo '$line' >> $new_file") if exists $namesH{$1}; 
	}
}
###copy result files ###
#blasts
my $new_blastDir = $new_dir."/results/blasts/$new_fam"; 
mkpath $new_blastDir;
my $blastFiles = $old_dir ."/results/blasts/$fam/*"; 
system("cp $blastFiles $new_blastDir/"); 
#clusters
my $clusterFiles = getAll::clustersFam($org, $class, $fam, 0);
if ($clusterFiles != 0){ 
	foreach my $cf (@$clusterFiles){
		my $full_cf = $old_dir ."/results/$cf";
		my $new_cf = $cf;
		my $old_taxa = $org."_".$class."_".$fam; 
		my $new_taxa = $new_org."_".$class."_".$new_fam; 
		$new_cf =~ s/$old_taxa/$new_taxa/; 
		my $new_full_cf = $new_dir."/results/$new_cf"; 
		open (my $cf_fh, $full_cf) or die "couldn't open $full_cf\n"; 
		open (my $new_cf_fh, ">>".$new_full_cf) or die "couldn't open $new_full_cf\n"; 
		my $cluster; 
		do{
			$cluster = readCluster($cf_fh); 
			(my $subfam) = $cluster =~ /[ACTG] name = \S+=(\S+)\n/; 
			print $new_cf_fh "$cluster" if $cluster ne "" and exists $namesH{$subfam}; 
		}while ($cluster ne ""); 
		close($new_cf_fh); 
		close($cf_fh); 
	}
}
else{
	print "no clusterFiles for $org $class $fam\n"; 
}

#read a cluster from an open file handle of a cluster file
sub readCluster{
	my $fh = shift; 
	my $line = '';
	my $cluster = ''; 
	#read cluster
	do{
		$line = <$fh>;
		$cluster .= $line; 
	}while ($line ne '' and $line !~ /^End cluster\n$/); 
	
	#get source, target and subfam from cluster
	if ($cluster =~ /\S+/){ #non empty cluster
		return $cluster; 
	}
	else {
		return ""; 
	}
}
