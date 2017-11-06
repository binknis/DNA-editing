#Function: Filters a cluster file. 
#Main usage: to create new track-files based on a filtered cluster file
#Note: "GA" mismatch will always have G as source and A as target, even if we're searching for A-to-I (see variables: $MORE_DIV_ELEMENT and $MAPPING_ELEMENT)
use strict; 
use Data::Dumper;
 use Cwd;
###Multiple filters will be applied: 
## 1. Cleanness filter on cluster stats
#List of filters (values specified are defaults): 
	# Min. # mismatches per cluster  5
	# Min. # mismatch in (one or more) clusters  10
	# Percent of mismatch and reverse (e.g. GA / (GA + AG) * 100)  60
	# Percent of mismatch and complement (e.g. GA / (GA + CT) * 100)  60
	# Percent of all mismatches in alignment  33.333
	# Minimal span of clusters (bps)  60
my $NUM_MM_CLUSTERED = 10; 
my $PC_MM_VS_REVERSE = 60; 
my $PC_MM_VS_COMPLEMENT = 60;
my $PC_OF_ALL_MM = 33.333; 
my $MIN_SPAN_CLUSTERS = 60; 

## 2. "Divergence from consensus" filter (A-more-div)
#	Select clusters between elements where A-seq is more diverged than G-seq
my $RMSK_FILE_HAS_BIN_COL = 0; 
my $MORE_DIV_ELEMENT = "target"; #"source"
## 3. "Mapped to consensus" filter (map to G more than A)
#	
my $NUC_MAPPED_TO_CONS = "pre"; #"post" #e.g. for GA pre = map G; post = map A; 
my $MAPPING_ELEMENT = ($MORE_DIV_ELEMENT eq "source" ? "target" : "source"); #'source' if GA given for G>A and 'target' if GA given for A>G; Should be opposite than more_div_element
my $MIN_MAPPED_PRE_VS_POST_FRAC = 0.5; #min value (non-inclusive) that will be accepted for 


### Generally used data structures
#Per nucleotides
my @nucs = ("A", "C", "G", "T"); 
my %nucInd = (); 
for (my $i=0; $i < scalar(@nucs); $i++){ 
	$nucInd{$nucs[$i]} = $i; 
}
#Per mismatch
my @mms = ("GA", "CT", "GC", "GT", "CA", "TA", "AG", "TC", "CG", "TG", "AC", "AT"); #mm order in mmCountStr in clusterfile
my %mmInd = (); #Index of each mm in mm-array 
for (my $i=0; $i < scalar(@mms); $i++){ 
	$mmInd{$mms[$i]} = $i; 
}
my %rev = (); #reverse of each mismatch GA -> AG
my %comp = (); #complement of each mismatch GA -> CT
foreach my $m (@mms){
	(my $nuc1, my $nuc2) = split('', $m); 
	$rev{$m} = reverse $m;
	my $comp = $m; 
	$comp =~ tr/ACTG/TGAC/; 
	$comp{$m} = $comp; 
}

### Get and parse commandline input
(my $mm, my $clusterFile, my $outDir, my $rmsk_file, my $mapToConsensus_file) = @ARGV; 
$mm = "GA" unless $mm; 
$clusterFile = "/home/alu/binknis/Data/FULGL/LTR/results/Tracks/tracks_FULGL_LTR_ERVL_1e-0_5/GA/clusters_FULGL_LTR_ERVL_1e-0_5.tab" unless $clusterFile; 
$rmsk_file = "/home/alu/binknis/Avian/rawdata/rmsk/FULGL.rmsk.out" unless $rmsk_file; 

(my $trackDir, my $fileSuff) = $clusterFile =~ /(\S+)\/clusters_(\S+)\.tab/; 
###TEMP ***
$outDir = "/home/alu/binknis/tempFiltering"; mkdir $outDir;
#### END TEMP

#set the mismatch (either arg or get from directory of clusterfile)
unless ($mm){
($mm) = $clusterFile =~ /\/([ACTG][ACTG])\//;
}
print $mm ."\n"; 
#get pre and post mismatches 
my $mmPre, my $mmPost; 
if($MAPPING_ELEMENT eq "source"){
	($mmPre, $mmPost) = split('', $mm); 
}
else{
	($mmPost, $mmPre) = split('', $mm); 
}
my $nuc; 
if($NUC_MAPPED_TO_CONS eq "pre"){
	$nuc = $mmPre; 
} 
else{
	$nuc = $mmPost; 
}
$mapToConsensus_file = $trackDir ."/". "nucListPerPair_".$nuc."_".$fileSuff.".txt" unless $mapToConsensus_file; 

### Get data for "Divergence from consensus" filter (A-more-div)
my %coordsToDiv = (); 
open(my $rmsk_fh, "<", $rmsk_file) or die "open $rmsk_file\n"; 
while(my $l = <$rmsk_fh>){
	chomp $l; 
	next if $l =~ /^(\#|\s*(SW|score))|^\s*$/; #skip comments, headers and empty lines
	$l =~ s/^\s+//;
	$l =~ s/\s+/\t/g; 
	my @fs = split(/\t/, $l); #fields: 
	shift @fs if $RMSK_FILE_HAS_BIN_COL; #discard bin column
	# print "@fs\n"; 
	(my $swScore, my $milliDiv, my $milliDel, my $milliIns, 
		my $genoName, my $genoStart, my $genoEnd, my $genoLeft, my $strand, 
		my $repName, my $repClass, my $repFamily, my $repStart, my $repEnd, my $repLeft, my $id) = @fs; #
	$genoStart--; #convert to 0-base like in my DB (rmsk output is 1-base)
	$strand =~ s/C/-/; 
	my $coords = $genoName .":". $genoStart ."-". $genoEnd . $strand; 
	
	$coordsToDiv{$coords} = $milliDiv; 
}
close($rmsk_fh); 

# print Dumper(\%coordsToDiv); exit; #***

### Get data for Consensus mapping filter (Most G map G)
##input file format: 
#normal line - 
#scaffold11427:19468-20535-|scaffold2879:1963-3034+      0       0       5       0       0       0       1       0
#line where none were successfully mapped - 
#scaffold12918:46444-46943+|scaffold45232:13026-13546-   0       0       0       0
my %pairsMapped = (); #pairs that passed the mapping filter
my $preInd = $nucInd{$mmPre}; 
my $postInd = $nucInd{$mmPost}; 
open(my $mapCons_fh, $mapToConsensus_file) or die "open $mapToConsensus_file\n"; 
while(my $l = <$mapCons_fh>){
	chomp $l; 
	my @fs = split(/\t/, $l); #fields: 
	next if scalar(@fs) < 9; #lines that weren't mapped should be skipped
	next if $fs[1 + $preInd]==0; #avoid division by 0
	#insert only pairs passing filter into hash
	if ( $fs[1 + $preInd] / ($fs[1 + $postInd] + $fs[1 + $preInd]) > $MIN_MAPPED_PRE_VS_POST_FRAC ){
		$pairsMapped{$fs[0]} = 1; 
	}
}
close($mapCons_fh);

print Dumper(\%pairsMapped); exit; #***

### Create and open output files
my $filter1 = "cleanAlign"; 
my $filter2 = $MORE_DIV_ELEMENT."MoreDiv"; 
my $filter3 = "moreMapPre";
my $clusters_filter1_outfile = $outDir ."/". "clusters_".$fileSuff."_".$filter1.".tab"; 
my $clusters_filter2_outfile = $outDir ."/". "clusters_".$fileSuff."_".$filter1."_".$filter2.".tab"; 
my $clusters_filter3_outfile = $outDir ."/". "clusters_".$fileSuff."_".$filter1."_".$filter2."_".$filter3.".tab"; 
open(my $f1_fh, ">".$clusters_filter1_outfile) or die "open $clusters_filter1_outfile\n"; 
open(my $f2_fh, ">".$clusters_filter2_outfile) or die "open $clusters_filter2_outfile\n"; 
open(my $f3_fh, ">".$clusters_filter3_outfile) or die "open $clusters_filter3_outfile\n"; 

### Read and filter cluster file
open(my $clusts_fh, $clusterFile) or die "open $clusterFile\n"; 
while(my $l = <$clusts_fh>){
	chomp $l; 
	(my $mismatch, my $assembly, my $class, my $family, my $subfam, 
		my $coordsS, my $coordsT, my $num_this_mm, my $num_all_mms, 
		my $total_prob, my $whereS, my $whereT, my $num_clusts, 
		my $mmSerials, my $clusters_span_woGaps, my $clusters_span, my $alignment_len, my $mmCount_str) = split(/\t/, $l); #fields: 
	
	#get num each mutation in alignment
	my @mmCounts = split('|', $mmCount_str);
	
	my $num_mm_clustered = scalar(split(',', $whereS)); 
	my $pc_mm_vs_reverse = $num_this_mm / ($num_this_mm + $mmCounts[$mmInd{$rev{$mismatch}}]) * 100; 
	my $pc_mm_vs_complement = $num_this_mm / ($num_this_mm + $mmCounts[$mmInd{$comp{$mismatch}}]) * 100; 
	my $pc_of_all_mm = $num_this_mm / $num_all_mms * 100; 
	
	#apply filter 1 - cleanness
	if($num_mm_clustered < $NUM_MM_CLUSTERED or 
		$pc_mm_vs_reverse < $PC_MM_VS_REVERSE or 
		$pc_mm_vs_complement < $PC_MM_VS_COMPLEMENT or 
		$pc_of_all_mm < $PC_OF_ALL_MM or 
		$clusters_span < $MIN_SPAN_CLUSTERS){ #skip lines that don't pass filter
		next; 
	}
	
	print $f1_fh $l ."\n"; #print clusters that past filter1
	
	#apply filter 2 - A more div
	if($MORE_DIV_ELEMENT eq "target"){
		# print "here1 " . $coordsToDiv{$coordsS} . " <-source, target-> " . $coordsToDiv{$coordsT} ."\n"; 
		# print $coordsS ."\t". $coordsT ."\n"; 
		if($coordsToDiv{$coordsS} > $coordsToDiv{$coordsT}){ #***check if need >=
			next; 
		}
	}
	else{
		# print "here2 " . $coordsToDiv{$coordsS} . " <-source, target-> " . $coordsToDiv{$coordsT} ."\n"; 
		# print $coordsS ."\t". $coordsT ."\n"; 
		if($coordsToDiv{$coordsS} < $coordsToDiv{$coordsT}){ #***check if need >=
			next; 
		}
	}
	print $f2_fh $l ."\n"; #print clusters that past filter1 and filter2
	
	#apply filter 3 - More map G
	next unless exists $pairsMapped{$coordsS ."|". $coordsT};
	print $f3_fh $l ."\n"; #print clusts passed all 3 filters
	
}
close($clusts_fh); 

close($f1_fh); 
close($f2_fh); 
close($f3_fh); 
	