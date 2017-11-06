package ProcessAlignmentByLength;
use strict;
use FindClustersByLength;

sub process_alignment_by_length
{
	(my $hsp_ref, my $query_name, my $hit_name, my $allMMs, my $cmd_line_ref) = @_; 
	my $hsp = $$hsp_ref;
	my $len = int( ( $hsp->length('query') + $hsp->length('hit') ) / 2 );

	if ( $hsp->strand('query') != $hsp->strand('hit') )
	{
		print
		"Alignement $query_name<=>$hit_name between different strands!!!\n";
		return;
	}
	
	my %no_rev = ("ga" => 0, "ct" => 0, "gc" => 0, "gt" => 0, "ca" => 0, "ta" => 0); 
	my %rev = ("ag" => 0, "tc" => 0, "cg" => 0, "tg" => 0, "ac" => 0, "at" => 0); #these should be reversed in output, because I want only one direction between each pair of nucs (ga, ct, gc, gt, ca, ta).
	
	my %transition = ("ag" => 1, "tc" => 1, "cg" => 0, "tg" => 0, "ac" => 0, "at" => 0, 
				"ga" => 1, "ct" => 1, "gc" => 0, "gt" => 0, "ca" => 0, "ta" => 0); #transitions and transversions are handeled differently after parsing
	
	my %mms = (); #position arrays for post-editing sequence (=target, child)
	my %mms_parent = (); #position arrays for pre-editing sequence (=source, parent)
	my %inds = (); #serial numbers of mismatch indices
	
	my %counts = (); 
	my %counts_rev = (); 
	foreach my $m (keys %transition){
		$counts{$m} = 0; 
		$counts_rev{$m} = 0; 
	}
	
	
	my $str1      = lc $hsp->query_string;
	my $str2      = lc $hsp->hit_string;
	my $alignment = $hsp->homology_string;

	my @s1 = split( //, $str1 );
	my @s2 = split( //, $str2 );
	
	my $q_start  = $hsp->start('query');
	my $s_start  = $hsp->start('subject');
	my $total_mms = 0;
	my $q_gaps   = 0;
	my $s_gaps   = 0;
	my $ag_count = 0;
	my $tc_count = 0;
	while ( $alignment =~ /( )/g )
	{
		my $p = ( pos $alignment ) - 1; #pos function is 1-base (counts from zeroth position). 
		my $m = $s1[$p] . $s2[$p]; #the mismatch
		if ( $s1[$p] eq '-' ) #query gap
		{
			$q_gaps++;
			next;
		}
		elsif ( $s2[$p] eq '-' ) #subject gap
		{
			$s_gaps++;
			next;
		}
		elsif (exists $no_rev{$m}){ #correct order (s1 is pre-edited nuc coined 'parent'; e.g. 'g' in g>a editing)
			push ( @{$mms{$m}}, $s_start + $p - $s_gaps ); #mms in 
			push ( @{$mms_parent{$m}}, $q_start + $p - $q_gaps ); #mms in parent
			push( @{$inds{$m}}, $total_mms);
			$counts{$m}++; 
		}
		elsif (exists $rev{$m}) { #reversed order (s2 is actual parent)
			push ( @{$mms{$m}}, $q_start + $p - $q_gaps ); #mms in 
			push ( @{$mms_parent{$m}}, $s_start + $p - $s_gaps ); #mms in parent
			push( @{$inds{$m}}, $total_mms );
			$counts{$m}++;
		}
		else{ #some non-standard nuc character (e.g. N) - do nothing (but total_mms will be incremented
		}
		$total_mms++;
	}	
	
	#reverse counts_rev hash - needed for command with reversed parent/edited 
	foreach my $m (keys %counts){
		$counts_rev{(reverse $m)} = $counts{$m}; 
	}
		
	## Calculate background probabilities for mismatches
	my %probs = (); 
	#calc transition probs
	$probs{"ct"} = $probs{"tc"} = ($counts{"ga"} + $counts{"ag"}) / ( 2 * $len ); #Each CT and TC probabilities are calculated by half of complementary mismatches
	$probs{"ga"} = $probs{"ag"} = ($counts{"ct"} + $counts{"tc"}) / ( 2 * $len ); #vice versa
	
	#calc transversion probs
	if ($allMMs){ 
		my $transition_sum = $counts{"ga"} + $counts{"ag"} + $counts{"ct"} + $counts{"tc"}; 
		foreach my $m (keys %transition){ #each type of mismatch (not only transitions...)
			next if $transition{$m}; #skip transitions
			$probs{$m} = ($total_mms - $transition_sum - $counts{$m} - $counts{(scalar reverse $m)}) / (6 * $len); #6 other transitions are basis for probability of this transition (all other transitions except for the rev of this one)
		}
	}	
	
	## Find clusters
	foreach my $m (keys %transition){ #iterates ALL types of mms
		if (not $allMMs){ #if only transitions - skip nontransitions
			next if $m !~ /^(ag|ga|ct|tc)$/;
		}
		if (exists $no_rev{$m}){ #correct order (e.g. ga)
			FindClustersByLength::find_clusters_by_length
					( \@{$inds{$m}}, \@{$mms{$m}}, \@{$mms_parent{$m}}, $query_name, $hit_name, $probs{$m}, $m, $total_mms, $len, $cmd_line_ref,  \%counts); #2 is flag for transversion
		}
		else{ #reversed (e.g. ag) -- #Note: no need to reverse mms and mms_parent because it is reversed above, upon alignment parsing
			FindClustersByLength::find_clusters_by_length
					( \@{$inds{$m}}, \@{$mms{$m}}, \@{$mms_parent{$m}}, $hit_name, $query_name, $probs{$m}, (scalar reverse $m), $total_mms, $len, $cmd_line_ref, \%counts_rev); #mm arrays and query/hit are reversed ; 2 is flag for transversion;
		}
	}
}

1;
