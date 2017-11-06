#Function: Download chain files for all possible pairs in a list of assemblies (for liftover)
#Input: 1. assemblies - file containing a list of assembly files (in one line or multiple lines) or a comma-delimited list of assemblies
#		2. outdir - output directory for chain files
use strict; 
use lib $ENV{HOME} . "/Perl_DNAE"; 
use getAll; 
use Cwd;
(my $assemblies, my $outdir) = @ARGV; 
$outdir = cwd() unless $outdir; 
#get list fo assemblies from input list or file
my @as = (); 
if($assemblies){ #comma-delim list
	@as = split(',', $assemblies); 
}
elsif (-e $assemblies){ #file name
	my $assemblyFile = $assemblies; 
	my $lines = getAll::lines($assemblyFile); 
	foreach my $l (@$lines){
		chomp $l; 
		push (@as, split(/\s+/, $l)); 
	}
}
else{
	die "bad input for assemblies\n"; 
}

#Download chain files for list of assemblies
foreach my $from (@as){
	foreach my $to (@as){
		next if $from eq $to; 
		(my $firstLetter, my $suffix) = $to =~ /(\S)(\S+)/; 
		my $ucAssembly = "" . (uc $firstLetter) . $suffix; 
		my $comm = "wget ftp://hgdownload.cse.ucsc.edu/goldenPath/$from/liftOver/".$from."To".$ucAssembly.".over.chain.gz"; 
		$comm =~ /liftOver\/(\S+)$/; 
		if ($outdir){
			my $outfile = "$outdir/$1";
			$comm .= " -O $outfile";
			my $outfileNoSuff = $outfile; 
			$outfileNoSuff =~ s/\.1$//; #there was an issue here that when downloaded into a dir containing the file then wget added ".1" at end of file name so the regulare "-e" may not work
			print $outfile ." ". $outfileNoSuff ."\n"; 
			next if(-e $outfile | -e $outfileNoSuff); #file was already downloaded
		}
		# print "not nexting\n"; 
		system($comm);
	}
}

