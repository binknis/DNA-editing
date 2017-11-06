#Function: Count amount of sequences per family in interval file exist after merging coords per family
#needed it for Nika project (To count non-redundant sequences). 
#Note: assumes interval file has 0-base starts
use strict; 

(my $interval_file) = @ARGV; 

(my $org) = $interval_file =~ /(\S+)\.interval/; 
my %lines_by_fam = (); 

open (INVL, $interval_file) || die "open $interval_file\n"; 
while (my $line = <INVL>){
	my @fields = split(/\s+/,$line); 
	$fields[6] =~ s/\?//; #erase question marks from family names (may never occur)
	$fields[5] =~ s/\?//; 
	$line =~ s/\?//g; 
	push(@{$lines_by_fam{$fields[5]."\t".$fields[6]}},$line); 
}
close(INVL); 

#Write output - file for each family
my $dir = "Dir_".$interval_file; 
mkdir $dir; 
foreach my $class_fam (sort keys %lines_by_fam){
	my $regex = "$class_fam".'\s*$'; 
	my $fam_file = $dir ."/".$class_fam.".interval"; 
	$fam_file =~ s/\t/_/; 
	open (FAM, ">$fam_file") || die "open $fam_file\n"; 
	print FAM @{$lines_by_fam{$class_fam}}; 
	close(FAM);
	system("sortBed -i $fam_file | mergeBed > $fam_file".".merged");
	open (MERGED, $fam_file.".merged") or die "open merged for $fam_file\n"; 
	my @lines = <MERGED>; 
	close(MERGED); 
	my $bps = getAmountBPs(\@lines);
	print $org ."\t". $class_fam ."\t". (scalar @lines)."\t". $bps."\n"; 
	unlink ($fam_file.".merged");
	unlink ($fam_file);
}
rmdir ($dir); 


### SUBROUTINES ###
sub getAmountBPs(){
	my $lines_ref = shift;
	my $bp_count; 
	foreach my $line (@$lines_ref){
		chomp $line; 
		my @fields = split(/\s+/,$line); 
		$bp_count += $fields[2] - $fields[1]; 
	}
	return $bp_count; 
}
