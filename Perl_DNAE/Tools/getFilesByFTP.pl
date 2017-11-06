#Function: Downloads files from UCSC for every subdirectory (e.g. assembly) listed
#you can enter "aaa" into a filename and it will be replaced with the assembly name you're trying to download
#file names containing the assembly name will be saved as is to pwd
#file names non containing the assembly name will be saved with a prefix "myAssembly_" 
#example for ftp path: "ftp://hgdownload.cse.ucsc.edu/goldenPath/gorGor3/database/ensGene.txt.gz"
use warnings;
use strict; 
use Net::FTP;

my @input = @ARGV; 

my $ftp_site = "hgdownload.cse.ucsc.edu"; 
#Note: assembly will usually be later inserted between the ftp site and the subdir. 
#my $subdir = "BigZips"; 
my $subdir = "database"; 
#my $file_to_download = "xenoRefGene.txt.gz"; 
#my $file_to_download = "ensGene.txt.gz"; 
my $file_to_download = "ensemblToGeneName.txt.gz"; 
my @fileList =(); 
if (@input){
	@fileList = @input; 
}
else{
	push(@fileList, $file_to_download); 
}

#my @assemblies = qw /proCap1/; 
my @assemblies = qw /ci2 strPur2 anoGam1 apiMel2 droAna2 droEre1 droGri1 dm3 droMoj2 droPer1 dp3 droSec1 droSim1 droVir2 droYak2 vicPac1 dasNov3 papAnu2 otoGar3 felCat4 panTro3 bosTau4 canFam2 turTru2 loxAfr3 nomLeu1 gorGor3 cavPor3 eriEur1 equCab2 hg19 dipOrd1 triMan1 calJac3 pteVam1 myoLuc2 micMur1 mm9 hetGla1 monDom5 ponAbe2 ailMel1 susScr2 ochPri2 ornAna1 oryCun2 rn4 rheMac2 proCap1 oviAri1 sorAra1 choHof1 saiBol1 speTri2 tarSyr1 sarHar1 echTel1 tupBel1 macEug2 cerSim1 cb3 ce6 priPac1 aplCal1 gadMor1 melUnd1 galGal3 latCha1 fr2 anoCar2 oryLat2 geoFor1 oreNil2 chrPic1 gasAcu1 melGal1 xenTro2 taeGut1 danRer7/; 
my @MissingAssemblies = (); 
my @ExistingAssemblies = (); 
my $ftp = Net::FTP->new($ftp_site);
$ftp->login("anonymous",'binknis@gmail.com');


foreach my $generic_file(@fileList){
	foreach my $assembly (@assemblies){
		my $dirToCheck = "/goldenPath/$assembly" ."/". $subdir;

		if ($ftp->cwd($dirToCheck))
		{
		  # print "Directory $dirToCheck exists\n"; 
		   push(@ExistingAssemblies, $assembly); 
		}
		else{
			#print "Directory $dirToCheck doesn't exist\n"; 
			push(@MissingAssemblies, $assembly); 
		}	
	}


	foreach my $assembly (@ExistingAssemblies){
		my $dir = "/goldenPath/$assembly" ."/". $subdir;
		my $final_generic = $generic_file; 
		$final_generic =~ s/aaa/$assembly/g; 
		my $saveTo = $assembly . "_". $final_generic unless $final_generic =~ /$assembly/; 
		$ftp->cwd($dir); 
		$ftp->get($generic_file, $saveTo) or print "error for $generic_file\n"; 
	}
}

$ftp->quit;