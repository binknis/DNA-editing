#Function: gets a sites gff file and removes redundancy
#Input: GFF file of sites created by createTrackFiles.pl 
#Output: 	1. GFF file, non-redundant sites. 
#			2. GFF file, redundant sites that were removed. 
#			3. GFF file, non-redundant sequences. 
#			4. GFF file, redundant sequences that were removed
#	PRINTS to stdout: The number of entries printed to each of the above files. (output format is self explanatory)
#Note: 	1. can accept header lines (will be skipped)
#		2. Can't accept gffs without +/-. 

( my $org, my $class, my $fam, my $pval, my $th, my $control) = @ARGV;
my $sites_gff; 
if ($org !~ /\.gff/){ #dir info input
	my $suffix = join ('_', $org,$class,$fam,$pval,$th).($control ? "_control" : ""); 
	$sites_gff = $ENV{"HOME"} ."/Data/$org/$class/results/Tracks/tracks_". $suffix ."/sites_A_".$suffix.".gff"; 
}
else{ #gff file input
	$sites_gff = $org; 
}
if (not -e $sites_gff){
	die "$sites_gff doesn't exist\n";
}


#(my $sites_gff) = @ARGV; 
#(my $sites_gff, my $min_sites_for_retaining_seq) = @ARGV; 
#$min_sites_for_retaining_seq = 0 if $min_sites_for_retaining_seq == 0; 
my %seqs = (); #coordinates of each sequence (original gff)
my %sites = (); #coordinates of each site (original gff) HoH - $sites{$site_coords}{$seq_coords}
my %seqs_nr = (); #non redundant sequences
my %seqs_r = (); #redundant seqs
my %sites_nr = (); #non redundant sites
my %sites_r = (); #redundant sites
my %lengths = (); #lengths of each sequece (for internal use)
my $prog_name = "DNAEfinder"; #const for GFF output
my $original_sites_count=0; 
#read input - original GFF
open(IN, $sites_gff) or die "open $sites_gff\n"; 
while (my $line = <IN>){
	chomp $line;
	my @fields = split (/\t/, $line);
	next if ($#fields != 8 || $fields[6] !~ /^[+-]$/); #skip header lines
	(my $chr, my $site, my $coords) = ($fields[0], $fields[3], $fields[8]);
	($chr, my $start, my $end, my $strand) = $coords =~ /([^:]+):(\d+)-(\d+)([+-])/;
	$seqs{$coords}=0;
	$lengths{$coords} = $end - $start + 1; 
	$sites{$chr.":".$site."-".$site.$strand}{$coords}=0;
	$original_sites_count++; 
}
close(IN); 

#associate each site with only one sequence (if site was initially redundant it's associated to the longest seq).
foreach my $site (keys %sites){
	my @containing_seqs = keys %{$sites{$site}}; 
	if (scalar @containing_seqs > 1){ #more than one sequence contains this site - select longest seq. 
		my $longest;
		my $max_len=0;
		foreach my $seq(@containing_seqs){ 
			if ($lengths{$seq} <= $max_len){ #save redundant site (option 1 - it wasn't even temporarily saved as longest)
				$sites_r{$site}{$seq}=0; 
			}
			if ($lengths{$seq} > $max_len){ 
				$sites_r{$site}{$longest}=0 if $max_len > 0; #save redundant site (option 2 - it was previously saved as longest)
				$longest = $seq;
				$max_len = $lengths{$seq}; 
			}
		}
		#save nr seq
		$sites_nr{$site}=$longest; 
		$seqs_nr{$longest}=0;
	}
	else{ #no redundancy for this site (contained in only 1 seq)
		$sites_nr{$site}=$containing_seqs[0]; 
		$seqs_nr{$containing_seqs[0]}=0; 
	}
}

### Print output ###

#print output - nr sites
my $sites_gff_nr = $sites_gff . ".nr"; 
my $nr_site_count = printGFF(\%sites_nr, $sites_gff_nr, \%sites_nr); 

#pring output - nr seqs
my $seqs_gff_nr = $sites_gff_nr;
$seqs_gff_nr =~ s/sites/seqs/; 
my $nr_seq_count = printGFF(\%seqs_nr, $seqs_gff_nr); 

##redundancy output

#print output - redundant sites
my $sites_gff_r = $sites_gff_nr;
$sites_gff_r =~ s/\.nr$/.redundant/;  
my $r_site_count = printGFF(\%sites_r, $sites_gff_r, \%sites_r); 

#output redundant sequences
foreach my $seq (keys %seqs){ #find redundant sequences (those that didn't have any sites uniquely associated with them)
	$seqs_r{$seq}=1 unless exists $seqs_nr{$seq};
}
my $seqs_gff_r = $sites_gff_r;
$seqs_gff_r =~ s/sites/seqs/; 
my $r_seq_count = printGFF(\%seqs_r, $seqs_gff_r); 

#print redundancy stats to STDOUT
print "Redundant in " . $sites_gff ."(sites and sequences):\t". $r_site_count ."\t". $r_seq_count ."\n"; 
print "Non-redundant in " . $sites_gff ."(sites and sequences):\t". $nr_site_count ."\t". $nr_seq_count ."\n"; 

#### SUBS ####

#Function: Print a set of coordinates to output
#$coords_ref: Hash or Array of coords for printing (if Hash then keys are coords)
#$group = info for 9th field (grouping) - 3 options: 1. hash to access each time with the key extracted from $coords_hash_ref; 2. HoH - a line of 1st key will be written for each 2nd key. 3. nothing (group will be count) 
#Note: the option to send an array as $coords_ref hasn't been invoked herein (but should work)
sub printGFF {
	(my $coords_ref, my $out_file, my $group) = @_; 
	open (COORDS, ">$out_file") || die "open $out_file\n";
	my $count=0;
	my $group_out;
	my @coords_ar; 
	#get coords from argument
	if ($coords_ref =~ /^ARRAY/){ #array of coords
		@coords_ar = sort @$coords_ref; 
	}
	elsif($coords_ref =~ /^HASH/){ #hash of coords
		@coords_ar = sort keys %$coords_ref; 
	}
	else{ #bad arg - shouldn't happen...
		print "bad input to printGFF\n"; 
	}
	#print output
	foreach my $coords(@coords_ar){
		my @groups;
		if($group =~ /^HASH/){ #group is a HASH ref
			if ($group->{$coords} =~ /^HASH/){ #HoH groups was arged - output line for each 2nd key
				foreach my $g (sort keys %{$group->{$coords}}){
					push (@groups, $g);
				}
			}
			else{ #regular Hash was arged - only 1 output line
				push (@groups, $group->{$coords});
			}
		}
		else{ #group wasn't arged - use simple count to seperate sequences from each other
			push (@groups, $count); 
		}
		(my $chr, my $start, my $end, my $strand) = $coords =~ /([^:]+):(\d+)-(\d+)([+-])/; #note: for sites, start == end.
		foreach my $group_out (@groups){ #only if $groups was HoH will this run more than once
			$count++; 
			print COORDS $chr."\t".$prog_name."\t"."."."\t".$start."\t".$end."\t"."."."\t".$strand."\t"."."."\t".$group_out."\n"; 
		}
	}
	close(COORDS); 
	return $count; #returns amount of lines printed
}