#If the name of the algorithm is changed from DNAEfinder then the grep needs to be modified as well.
#subfamFlag = if to create files for subfamilies
#Note make sure to create track files before running this! Otherwise the num of edited sites will be 0.


use strict; 
use getAll; 
# $GA = G or A 
( my $org, my $class, my $family, my $pval, my $th, my $GA, my $control, my $subdir) = @ARGV;

$subdir = 0 if $subdir eq ''; 
$control = 0 if $control eq '';
my $suffix = $org . "_" . $class . "_" . $family . "_" . $pval . "_" . $th; 
$suffix .= "_control" if ($control);

my $dir = "../Data/" . $org . "/" . $class . "/results"; 
$dir .= "/$subdir" unless ($subdir eq 0 or $subdir eq ""); 
$dir .= "/Tracks"; 
$dir .= "/tracks_" . $suffix;
my $siteFile = $dir ."/sites_". $GA ."_". $suffix .".gff";
my $siteListFile = $dir ."/siteList_". $GA ."_". $suffix .".txt";

print $org . "\t" . $class . "\t" . $family ."\t";
unless (-e $dir){ #print 0 and exit if trackFiles don't exist yet
	print "0\n"; 
	exit; 
}


system ("grep -c 'DNAEfinder' $siteFile");
#find by subfamily
my %sites = (); 
open( my $sites_handle, "<". $siteListFile ) or die "siteList file didn't open for $org $class $family: ";
while (my $line = <$sites_handle>)
{ 
	(my $name)  = ( $line =~ /^(\S+)\t/ );
	(my @where) = ( $'    =~ /(\d+)/g );
	(my $subfam) = $name =~ /\S+=(\S+)/;
	$sites{$subfam} += $#where + 1;
}
close($sites_handle);

#print subfam output
my $subfamFile = $dir ."/sf_siteCount_". $GA ."_". $suffix .".txt";
open (SUBFAM_SITES, ">".$subfamFile) || die "couldn't open $subfamFile\n"; 
foreach my $sf (sort keys %sites){
	print SUBFAM_SITES $org ."\t". $class ."\t". $family ."\t". $sf ."\t". $sites{$sf}."\n"; 
}
close(SUBFAM_SITES);
