use strict; 
use getAll;
use Math::BigInt; 
package getParams; 
# $clust_stats_file: full path.
#clust_stats file format: 
#Fields: Organism (0) Class (1)	Family(2)	Pval(3)	Threshold(4)	#editing(5)	#editing_control(6)	%FP_editing(7)	#edited(8)	#edited_control(9)	%FP_edited(10)
#($fam, $pval, $th, $editing, $editing_ctrl, $fp_editing, $edited, $edited_ctrl, $fp_edited)

#reads a cluster_stats file and finds the parameters which have a value in a column closest to "val"
#returns: 1-2. The closest set of parameters (pval, th); 3. the difference from "val".
#Note: "fields" are 0-based.
sub closestTo{
	(my $clust_stats_file, my $field, my $val) = @_;
	my $best_pval = -1; 
	my $best_th = -1; 
	my $best_diff = Math::BigInt->binf(); #pseudo infinity
	my $lines_ref = getAll::lines($clust_stats_file); 
	die "didn't open $clust_stats_file\n" unless $lines_ref;
	foreach my $line (@$lines_ref){
		 my @fields = split(/\s+/,$line); 
		if (abs ($fields[$field] - $val) < $best_diff){
			$best_diff = abs ($fields[$field] - $val); 
			($best_pval, $best_th) = $line =~ /\t1e-(\d+)\t(\d+)/;
		}
	}
	return ($best_pval, $best_th, $best_diff); 
}

#same as "closestTo" but returns the original line and not the parameters. 
sub closestTo_returnLine{
	(my $clust_stats_file, my $field, my $val) = @_; 
	my $best_line = ""; 
	my $best_diff = Math::BigInt->binf(); #pseudo infinity
	my $lines_ref = getAll::lines($clust_stats_file); 
	die "didn't open $clust_stats_file\n" unless $lines_ref;
	foreach my $line (@$lines_ref){
		 my @fields = split(/\s+/,$line); 
		if (abs ($fields[$field] - $val) < $best_diff){
			$best_diff = abs ($fields[$field] - $val); 
			$best_line = $line;
		}
	}
	return $best_line; 
}

sub closestToMaximizeBy{
	(my $clust_stats_file, my $field_close, my $val, my $field_max) = @_;
	my $best_pval = -1; 
	my $max_val = Math::BigInt->binf('-'); 
	my $best_th = -1; 
	my $best_diff = Math::BigInt->binf(); #pseudo infinity
	my $lines_ref = getAll::lines($clust_stats_file); 
	die "didn't open $clust_stats_file\n" unless $lines_ref;
	foreach my $line (@$lines_ref){
		 my @fields = split(/\s+/,$line); 
		if (abs ($fields[$field_close] - $val) < $best_diff){
			$best_diff = abs ($fields[$field_close] - $val); 
			$max_val = $fields[$field_max];
			($best_pval, $best_th) = $line =~ /\t1e-(\d+)\t(\d+)/;
		}
		elsif (abs ($fields[$field_close] - $val) == $best_diff and $fields[$field_max] > $max_val){
			$max_val = $fields[$field_max]; 
			($best_pval, $best_th) = $line =~ /\t1e-(\d+)\t(\d+)/;
		}
	}
	return ($best_pval, $best_th, $best_diff);
}

#This function's if-elses are different than the above. 
#This is the one for preferred use.
#Finds: The params within margin (target +/- margin) from target value, which maximize the "field_max" column. 
#Note: if no params within margin were found - returns "Family\tNA"

sub closestToMaximizeBy_returnLine{
	(my $clust_stats_file, my $field_close, my $val, my $field_max, my $diff_margin) = @_;
	my $best_line = ""; 
	my $max_val = Math::BigInt->binf('-');
	my $best_diff = Math::BigInt->binf();
	my $lines_ref = getAll::lines($clust_stats_file);
	die "Didn't open $clust_stats_file\n" unless $lines_ref;
	foreach my $line (@$lines_ref){
		my @fields = split(/\s+/,$line); 
		my $diff = abs ($fields[$field_close] - $val);
		my $val = $fields[$field_max]; 
		next if $diff > $diff_margin; 
		if ($diff < $best_diff || ($diff == $best_diff and $val > $max_val) ){ #best_diff is out of margin and a smaller diff is found
			$best_diff = $diff;
			$max_val = $val; 
			$best_line = $line;
		} 
		# print $line; #***
	}
	#return best line or NA if none were within margin. 
	if ($best_line eq ""){
		(my $fam) = $lines_ref->[0] =~ /^(\S+)/; 
		return $fam ."\t". "NA" . "\n";
	}
	return ($best_line);
}

#returns hash with keys as "$pval=$th", if constraint was met. 
sub nonZero{
	(my $clust_stats_file, my @field_nums) = @_; 
	my %nonZeroParams = (); 
	my $zero_found; 
	my $lines_ref = getAll::lines($clust_stats_file); 
	die "didn't open $clust_stats_file\n" unless $lines_ref;
	foreach my $line (@$lines_ref){
		$zero_found = 0; 
		my @fields = split(/\s+/,$line); 
		foreach my $val (@field_nums){
			$zero_found = 1 if $fields[$val] == 0; 	
		}
		
		unless ($zero_found){
			(my $pval, my $th) = $line =~ /\t1e-(\d+)\t(\d+)/;
			$nonZeroParams{"$pval=$th"}=1; 
		}
	}
	return \%nonZeroParams; 
}


#This function's if-elses are different than the above. 
#This is the one for preferred use.
#Finds: The params within margin (target +/- margin) from target value, which maximize the "field_max" column. 
#Note: if no params within margin were found - returns "Family\tNA"

sub inMarginMaximizeBy_returnLine{
	(my $clust_stats_file, my $field_close, my $val, my $field_max, my $diff_margin) = @_;
	my $best_line = ""; 
	my $max_val = Math::BigInt->binf('-');
	my $best_diff = Math::BigInt->binf();
	my $lines_ref = getAll::lines($clust_stats_file);
	die "Didn't open $clust_stats_file\n" unless $lines_ref;
	foreach my $line (@$lines_ref){
		my @fields = split(/\s+/,$line); 
		my $diff = abs ($fields[$field_close] - $val);
		my $val = $fields[$field_max]; 
		if ($diff <= $diff_margin and $val > $max_val){ #best_diff is out of margin and a smaller diff is found
			$best_diff = $diff;
			$max_val = $val; 
			$best_line = $line;
		} 
		# print $line; #***
	}
	#return best line or NA if none were within margin. 
	if ($best_line eq ""){
		(my $fam) = $lines_ref->[0] =~ /^(\S+)/; 
		return $fam ."\t". "NA" . "\n";
	}
	return ($best_line);
}
















return 1; 