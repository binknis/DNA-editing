#Function: Package of subroutine "find_clusters_by_length" which uses 2 subs: "max" and "write_cluster".  
#find_clusters_by_length: 
#max: 
#write_cluster:

#Variables: 
#th := threshold.

#Notes after changes for screen: 
#1. add pval_l/h and th_l/h
#2.  pval_h: loose (1e-6). pval_l: stringent (1e-15).
#3. example of th values: th_l == 5. th_h == 12. 

#Notes for future optimization: 
# possible to use refs to arrays instead of copying them

#Output Notes: 
#	1. write_cluster_tabular writes tabular output. writes also GA and CT mismatches to ga and ct subdirs, respectively; If you want to go back to older version, replace the function call to write_cluster (see changes below)

package FindClustersByLength;
use strict;
use Math::NumberCruncher; 
my $tabular; 
my @mm;
my @mm_parent; 
my @ind;
my %whereHash;
my $num_all_mms;
my $alignment_len;
my $organism; my $class; my $family;
my $pval_l; my $pval_h; 
my $th_l; my $th_h; 
my $s_name; my $t_name; #source and target (in G>A editing, G is source and A is target)

my @mmTypes = qw/ga ct gc gt ca ta ag tc cg tg ac at/; 

sub find_clusters_by_length
{	
	( my $ind_ref, my $mm_ref, my $mm_parent_ref, $s_name, $t_name, my $prob, my $mmType, $num_all_mms, $alignment_len, my $cmd_line_ref, my $mmCounts)	= @_; 
	@ind = @$ind_ref; #serial mm number
	@mm = @$mm_ref;  #pos in child (A) for G>A editing
	@mm_parent = @$mm_parent_ref; #pos in parent (G) for G>A editing
	
	( $organism, $class, $family, $pval_h, $pval_l, $th_l, $th_h , $tabular ) = @$cmd_line_ref; #pval and th order example: 5 16 5 16 (latter are 1e-5 to 1e-16)

	return if ($#mm + 1 < $th_l); #less mismatches than threshold - no possible clusters!
	my %num_clusts = (); 
	my %total_prob = (); 
	my @pvals = (); 
	
	#create pval array (loose to stringent)
	for (my $i=$pval_h; $i<=$pval_l ;$i++){
		push (@pvals, "1e-$i"); 
	}
	
	#create params' variables
	foreach my $pval (@pvals){
		for (my $th=$th_l; $th<=$th_h; $th++){
			$total_prob{$pval}{$th} = 1;
			$num_clusts{$pval}{$th} = 0; 
			$whereHash{$pval}{$th} = (); 
		}
	}
	
	#Get START and END indices of all consecutive G->A mismatch intervals
	my $start; 
	my $end; 
	my @starts; 
	my @ends; 
	my $i=0; 
	my $j;
	while ($i < $#ind-$th_l+2){ #maximal window left in mm array is not smaller than th_l
		$start = $i; 
		while ($i < $#ind && $ind[$i]+1 == $ind[$i+1]){ #increment index past all consecutive G-to-A mismatches
			$i++;
		}
		$end = $i;
		if ($end-$start+1 >= $th_l){ #if interval is sufficiently long (at least for th_l) - add the interval's start & end
			push(@starts, $start);
			push(@ends, $end);
		}
		$i++;
	}
	return if $#starts == -1; #no sufficiently long intervals found
	
	### SEARCH INTERVALS FOR CLUSTERS ###
	my $size; my $max_th;
	my $clust_count = 0; 
	my $old_clust_count;
	for (my $interval=0; $interval<=$#starts; $interval++){ #for each interval to search in
		foreach my $pval (@pvals){ #each pval (loose to stringent)
			$old_clust_count = $clust_count;
			$start = $starts[$interval];
			$end = $ends[$interval];
			$i = $start;
			$size = $end - $start + 1;
			until ($size < $th_l){ 
				$max_th = min($th_h, $size);
				$j = $i + $size - 1;
				my $cprob=0;
				if ($pval < 1){ #any pval that isn't 1 (1e-0; used to accept any cluster disregarding pval)
					$cprob = Math::NumberCruncher::Binomial( $mm[$j] - $mm[$i] + 1, $j - $i + 1, $prob );
				}
				if ($cprob < $pval){
					for(my $th=$th_l; $th<=$max_th; $th++){ #***write by pval?
						push( @{$whereHash{$pval}{$th}}, ($i..$j) );
						$total_prob{$pval}{$th} *= $cprob; 
						$num_clusts{$pval}{$th}++;
						$clust_count++;
					}
					#found cluster, thus finished searching for this pval 
					last; #proceed to next pval
				}
				else{ #advance beginning of window unless window reached end, in which case window-size is decremented, and window is shifted back to start
					if($j < $end){
						$i++;
					}
					else{
						$size--;
						$i=$start;
					}
				}
			}
			if ($clust_count - $old_clust_count == 0){ #no clusters found for permissant pval - no need to check stringent pvals.
				last; 
			}
		}
	}

	#write clusters
	return if $clust_count == 0; #no clusters found
	foreach my $pval(@pvals){
		for (my $th=$th_l; $th<=$th_h; $th++){
			if ($num_clusts{$pval}{$th} > 0){ #cluster found
				if ($tabular){ #tabular output. All files in mm-specific subdirs.
					write_cluster_tabular($pval, $th, $total_prob{$pval}{$th}, $num_clusts{$pval}{$th}, $mmType, $mmCounts);
				}
				else{ #Shai's original format. Only transversions are in subdirs.
					write_cluster($pval, $th, $total_prob{$pval}{$th}, $num_clusts{$pval}{$th}, $mmType);
				}
			}
		}
	}
}

#arg: array[2]. returns the highest value
sub max
{
	if ( $_[0] < $_[1] ) { return $_[1] }
	else { return $_[0] }
}

sub min
{
	if ( $_[0] > $_[1] ) { return $_[1] }
	else { return $_[0] }
}

#write clusters to a new file (Shai's original format)
#"_control" will be added to file name if it's a control
sub write_cluster
{	
	(my $pval, my $th, my $total_prob, my $num_clusts, my $mmType) = @_; 
	(my $preEdit, my $postEdit) = split('', uc $mmType); 
	my @where = @{$whereHash{$pval}{$th}};
	
	my $cluster_name;
	#create and open file
	if ($mmType eq "ga" or $mmType eq "ag") { #GA or AG (typical G>A DNA- or A>I RNA-editing) #Note: "ag" isn't supposed to happen because it's reversed to 'ga' before calling this command (same for "tc")
		$cluster_name = ">>../Data/$organism/$class/results/clusters_" . $organism . "_" . $class . "_"	. $family . "_"	. $pval . "_"	. $th . ".txt";
	}
	elsif ($mmType eq "ct" or $mmType eq "tc"){ #CT or TC (historically "control" for G>A DNA editing)
		$cluster_name =	">>../Data/$organism/$class/results/clusters_" . $organism . "_" . $class . "_" . $family . "_"	. $pval . "_"	. $th . "_control.txt";
	}
	else { #transversions (used to see background mutation rate)
		$cluster_name = ">>../Data/$organism/$class/results/$mmType/clusters_" . $organism . "_" . $class . "_"	. $family . "_"	. $pval . "_"	. $th . ".txt";
	}
	
	open( CLUSTS, $cluster_name );
	print CLUSTS "Found cluster:\n";
	print CLUSTS $preEdit . " name = $s_name\n";
	print CLUSTS $postEdit . " name = $t_name\n";
	print CLUSTS "Total probability = $total_prob, $num_clusts cluster";
	if ($num_clusts > 1) {print CLUSTS "s\n";} else {print CLUSTS "\n";}	
	my $str_format = '';
	foreach my $x (@where) { $str_format .= "%5d "; }
	$str_format .= "\n";
	print CLUSTS "Edited MM serial no.    : ";
	my @where1base;
	foreach my $i ( 0 .. $#where ) { $where1base[$i] = $where[$i] + 1; }
	printf CLUSTS $str_format , @where1base;
	print CLUSTS "Locations $preEdit seq: ";
	printf CLUSTS $str_format , @mm_parent[@where];
	print CLUSTS "Locations $postEdit seq: ";
	printf CLUSTS $str_format , @mm[@where];
	
	print CLUSTS "Total length: $alignment_len, Direct mismatches: ", scalar(@mm),	", All mismatches: $num_all_mms\n";
	print CLUSTS "End cluster\n";
	close(CLUSTS);
}

#Changes from "write_cluster": 
#	1. Writes also GA and CT mismatches to a respective subdir in "results"
#	2. All output is in tabular format - one tuple per alignment (HSP)
sub write_cluster_tabular
{	
	(my $pval, my $th, my $total_prob, my $num_clusts, my $mmType, my $mmCounts) = @_; 
		
	#Build tuple for output
	#The fields: assembly[1], taxa (org, class, fam, subfam)[2-5], mismatch-type, 
	(my $assembly, my $coordsS, $class, $family, my $subfam) = $s_name =~ /^\d+=([^=]+)=([^:]+:\d+-\d+[+-])=([^=]+)=([^=]+)=(\S+)/; 
	(my $coordsT) = $t_name =~ /^\d+=[^=]+=([^:]+:\d+-\d+[+-])=/; 

	my @where = @{$whereHash{$pval}{$th}};
	my @where1base;
	foreach my $i ( 0 .. $#where ) { $where1base[$i] = $where[$i] + 1; }
	my $mmSerials = join (',', @where1base); 
	my $whereS = join (',', @mm_parent[@where]); 
	my $whereT = join (',', @mm[@where]); 
	my $num_edited_sites = scalar(@mm); 
	
	my $clusters_span = $mm[$where[$#where]] - $mm[$where[0]] + 1; #num bps of region affected by editing (bps span from 1st mm of 1st cluster to last mm of last cluster)
	
	my $clusters_span_woGaps=0; 
	my $c_size; 
	my $pres=-1; #present
	for (my $i=1; $i<=$num_clusts; $i++){
		$pres++; 
		my $start = $pres;
		while ($pres < $#where and $where[$pres+1] == $where[$pres] + 1) { #next is successive serial number and haven't reached end of array
			$pres++; 
		};
		$c_size = $mm[$where[$pres]] - $mm[$where[$start]] + 1;
		$clusters_span_woGaps += $c_size; 
	}
	
	#str of count of each type of mismatch in alignment
	my $mmCount_str = ""; 
	foreach my $m (@mmTypes){
		$mmCount_str .= $mmCounts->{$m} . "|";
	}
	chop $mmCount_str; #erase trailing comma
	#create actual tuple for output
	my $tuple = (uc $mmType)."\t"
					.$assembly."\t".$class."\t".$family."\t".$subfam."\t"
					.$coordsS."\t".$coordsT."\t"
					.$num_edited_sites."\t".$num_all_mms."\t".$total_prob."\t"
					.$whereS."\t".$whereT."\t".$num_clusts."\t".$mmSerials."\t"
					.$clusters_span_woGaps."\t".$clusters_span."\t".$alignment_len."\t"
					.$mmCount_str;
	
	#create & open file & print output
	my $cluster_name;
	$cluster_name = ">>../Data/$organism/$class/results/". (uc $mmType) . "/clusters_" . $organism . "_" . $class . "_"	. $family . "_"	. $pval . "_"	. $th . ".tab";
	open( CLUSTS, $cluster_name );
	print CLUSTS $tuple ."\n"; 
	close(CLUSTS);
}

1;