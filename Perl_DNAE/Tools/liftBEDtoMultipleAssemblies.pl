use strict; 
use lib $ENV{HOME} . "/Perl_DNAE"; 
use getAll; 
use Data::Dumper;

#CONSTS
my $noOverrideLift = 1; #if the liftoverfile exists - don't run the command again

(my $bedAssembly, my $bedFile, my $assemblies, my $outDir, my $outSuffix, my $chainDir) = @ARGV; 
$chainDir = "/home/alu/binknis/binknis_data/DBs/UCSC/liftOver/over.chains" unless $chainDir; 

#get assemblies to lift to (NOTICE: assems and numAssems contains target assemblies only!).
(my $assems_r, my $numAssems ) = getAssemblies($assemblies, $bedAssembly); 
my @assems = @$assems_r; 

#lift all sites in input BED file to other assemblies
liftToMultipleAssemblies($bedAssembly, $assems_r, $chainDir, $outDir, $outSuffix, $noOverrideLift); 

###convert all results into table
#get all results into hash
my %liftRes = (); #$liftCoord{$toAssembly}{NA|0|1|2} #NA is initialization (pipes resemble OR)
my %liftCoord = (); #$liftCoord{$toAssembly}{NA|coords} #NA is initialization (pipe resembles OR)
#init array in hash for each coord in input file
initLiftArray($bedFile, \%liftRes, \%liftCoord, $numAssems); 


#get liftover results (lift, no-lift, partial)
getLiftRes($bedAssembly, $assems_r, $outDir, $outSuffix, \%liftRes, \%liftCoord); 

printLiftTable($bedAssembly, \@assems, $outDir, \%liftCoord, \%liftRes);
# printLiftTable($bedAssembly, \@assems, $outDir, \%liftCoord, \%liftRes, "tempLiftTable"); #TEMP for comparison to original.

###Function: get list of assemblies to lift to
#Input: $assemblies - either a comma-delimited list of assemblies or a file with a list of assemblies (whitespace delimited; can be multiline).
#Output: 	1. ref to array of list of assemblies
#			2. number of assemblies
sub getAssemblies{
	(my $assemblies, my $assemToRemove) = @_; 
	my @assems = (); 
	if($assemblies =~ /,/){ #comma-delim list
		@assems = split(',', $assemblies); 
	}
	elsif (-e $assemblies){ #file name
		my $assemblyFile = $assemblies; 
		my $lines = getAll::lines($assemblyFile); 
		foreach my $l (@$lines){
			chomp $l; 
			push (@assems, split(/\s+/, $l)); 
		}
	}
	else{
		die "bad input for assemblies\n"; 
	}
	
	@assems = grep { $_ ne $assemToRemove } @assems if($assemToRemove); #remove bedAssembly from assembly list

	my $numAssems = scalar(@assems);
	return (\@assems, $numAssems)
}

#Function: run liftover for the bedAssembly file to multiple assemblies
sub liftToMultipleAssemblies{
	(my $bedAssembly, my $assems_ref, my $chainDir, my $outDir, my $outSuffix, my $noOverrideLift) = @_; 
	foreach my $as (@$assems_ref){
		next if $as eq $bedAssembly; #don't lift to self...
		(my $firstLetter, my $suffix) = $as =~ /(\S)(\S+)/; 
		my $fromTo = $bedAssembly . "To" . (uc $firstLetter) . $suffix;
		my $chain =  $chainDir ."/". $fromTo . ".over.chain.gz";
		next unless -e $chain; #skip if chain file doesn't exist
		my $out = $outDir ."/". $fromTo . ($outSuffix ? ".". $outSuffix : "") . ".lifted"; #/home/alu/binknis/RNA_med_DNAE/Data/Alu/liftRes
		my $unlifted = $outDir ."/". $fromTo . ($outSuffix ? ".". $outSuffix : "") . ".unlifted"; 
		my $comm = "/home/alu/binknis/binknis_data/DBs/UCSC/liftOver/liftOver $bedFile $chain $out $unlifted"; 
		if(-e $out and $noOverrideLift){
			next; 
		}
		elsif(-e $chain){
			# print $comm ."\n"; 
			system ($comm); # liftOver input.bed hg18ToHg19.over.chain.gz output.bed unlifted.bed
		} 
		else{
			print $chain . " doesn't exist\n"; 
		}
	}
}

#Function: init lift result arrays
sub initLiftArray{
	(my $bedFile, my $liftRes_r, my $liftCoord_r, my $numAssems) = @_; 
	my $lines = getAll::lines($bedFile); 
	foreach my $l(@$lines){
		next unless $l =~ /\S+/; #skip empty lines in file
		chomp $l; 
		my @f = split("\t", $l); 
		$liftRes_r->{$f[3]} = (); 
		$liftCoord_r->{$f[3]} = (); 
		for(my $i=0; $i<$numAssems; $i++){
			push(@{$liftRes_r->{$f[3]}}, "NA"); #0-deleted, 1-successful liftover, 2-partial liftover or NA (init value)
			push(@{$liftCoord_r->{$f[3]}}, "NA"); #coordinates or NA (init value)
		}
	}
}

#Function: insert liftover results to result Hashes (liftRes and liftCoord)
sub getLiftRes{
	(my $bedAssembly, my $assems_r, my $outDir, my $outSuffix, my $liftRes_r, my $liftCoord_r) = @_; 
	my $asCount = -1; #for accession if index in hash of array
	for my $as (@$assems_r){
		$asCount++; 
		(my $firstLetter, my $suffix) = $as =~ /(\S)(\S+)/;
		my $fromTo = $bedAssembly . "To" . (uc $firstLetter) . $suffix;
		
		#get lifted coords
		my $out_lifted = $outDir ."/". $fromTo . ($outSuffix ? ".". $outSuffix : "") . ".lifted";
		open(my $l_fh, $out_lifted) or die "open $out_lifted\n"; 
		while(my $l = <$l_fh>){
			next unless $l =~ /\S+/; #skip empty lines in file
			chomp $l; 
			my @f = split("\t", $l); 
			my $from = $f[3]; 
			my $to = $f[0] .":". $f[1] ."-". $f[2] . $f[5]; 
			$liftRes_r->{$from}[$asCount] = 1; #value for successfully lifted
			$liftCoord_r->{$from}[$asCount] = $to; 
		}
		close($l_fh); 
		
		#get partially lifted
		my $out_unlifted = $outDir ."/". $fromTo . ($outSuffix ? ".". $outSuffix : "") . ".unlifted"; 
		open(my $ul_fh, $out_unlifted) or die "open $out_unlifted\n"; 
		while(my $l = <$ul_fh>){
			next unless $l =~ /\S+/; #skip empty lines in file
			chomp $l; 
			if($l =~ /^#/){
				my $liftVal = ($l =~ /^#Deleted/ ? 0 : 2); #0- flag for deleted, 2- flag for partially deleted
				$l = <$ul_fh>; #get line with coords 
				chomp $l; 
				my @f = split("\t", $l); 
				my $from = $f[3]; #from is in 4th col in this file
				$liftRes_r->{$from}[$asCount] = $liftVal; #value for successfully lifted
			}
		}
		close($ul_fh); 
	}
}

#Function: Print output: print liftTable to file
sub printLiftTable{
	(my $bedAssembly, my $assems_r, my $outDir, my $liftCoord_r, my $liftRes_r, my $out_prefix) = @_; 
	$out_prefix = "liftTable" unless $out_prefix; #set default output
	my $outTable = $outDir ."/".$out_prefix."_".$bedAssembly."_To_".join('_', @$assems_r).".tab";
	open(my $tab_fh, ">". $outTable) or die "open $outTable\n"; 
	print $tab_fh join("\t", "coords", (map { "res_" . $_ } @$assems_r),  (map { "coords_" . $_ } @$assems_r) ) ."\n"; 
	foreach my $c (sort keys %$liftCoord_r){
		next unless $c =~ /\S+/; #skip empty lines
		print $tab_fh $c ."\t".  join("\t", @{$liftRes_r->{$c}}) ."\t". join("\t", @{$liftCoord_r->{$c}}) . "\n";
	}
	close($tab_fh);  
}
