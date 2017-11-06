#Function: creates a fasta file from a list of requested sequences from the DB (in Data dir)
#Input: coordinate file (can be multiple format, but must have coords as 'chr:start-end' and afterwards 'subfamily')
#		column numbering is 0 based 
# $fastaOut : can provide with or without suffix (.fa will be the suffix; and _subfam.fa if file per subfam is activated)
#			  Can contain full path!
# $outFormatFlag : 0 - full, 1 - only coords, 2 - "coords=subfam", HASH - coords-to-new-defline hash
# $filePerSubFamOutput : 0 - one large file. 1 - a file per subfamily
# $columnsToGet : specifies the column containing the coords; 0-base.
# $filePerFamInput : 1 - DB contains single files per. 0 - DB contains file per subfam. 
#Output: Fasta file with sequences (1 per family if flag is on)
#Notes: 1. coord filename must contain the suffix created by createTrackFiles.pl (Org, class, fam, pval, th)
#		2. subfam names with slashes and/or question-marks will have them retained in deflines but are substituted with underscores in filenames. 
use strict; 
use File::Path qw(mkpath);
use Bio::SeqIO;
(my $coordFile, my $fastaOut, my $columnToGet, my $outFormatFlag, my $filePerSubFamOutput, my $coordsToSubfam) = @ARGV; 
my $coords; 
my $coordToGet; 
my $org, my $class, my $fam; 
my %subfamToCoords = (); 
#extract Organism, Class and family from file name
if ($coordFile =~ /\S+[\/]*[^_\/]+_([^_\/]+)_([^_\/]+)_(\S+)_(1e-\d+)_(\d+)(_control)?.txt$/){
	($org, $class, $fam) = ($coordFile =~ /\S+[\/]*[^_\/]+_([^_\/]+)_([^_\/]+)_(\S+)_(1e-\d+)_(\d+)(_control)?.txt$/); 
}
else{
	die "insufficient data in file name: $coordFile - need Organism, Class and Family\n"; 
}	

#create hash with the subfamily and coords wanted ({subfam}{coords})
open (COORDS, "$coordFile") or die "$coordFile didn't open\n"; 
while (my $line = <COORDS>){
	chomp $line; 
	my @fields = split (/\s+/,$line); 
	$coordToGet = $fields[$columnToGet];
	($coords, my $subfam) = $coordToGet =~ /([^=]+:\d+-\d+[+-])\S*=(\S+)$/;
	if($coordsToSubfam){ #replace subfam found in defline with subfam in $coordsToSubfam HASH
		$subfam = $coordsToSubfam->{$coords}; 
	}
	$subfam =~ s/\?//g;  
	$subfam =~ s/\//_/g;
	$subfamToCoords{$subfam}{$coords} = 0;
}
close (COORDS);

#copy the wanted sequences from the subfamily files in the DB
my $outStream;
$fastaOut =~ s/\.fa$//; #to avoid adding .fa to existent .fa suffix
$outStream = Bio::SeqIO->new( -file => ">".$fastaOut.".fa",  -format => 'fasta' ) unless $filePerSubFamOutput;
foreach my $subfam (sort keys %subfamToCoords){
	my $dbFile = "../Data/".$org."/".$class."/db/files_".$fam."/Seq_".$subfam; 
	if($coordsToSubfam){ #only one file per family in input DB
		"../Data/".$org."/".$class."/db/files_".$fam."/Seq_".$fam; 
	}
	if ($filePerSubFamOutput){
		$outStream = Bio::SeqIO->new( -file => ">".$fastaOut."_".$subfam .".fa",  -format => 'fasta' ); #subfamily-specific output file
	}
	copySubfamSequences($dbFile, $outStream, $subfamToCoords{$subfam}, $outFormatFlag, $coordsToSubfam); 
}

#Function: reads the fasta db file and copies only those existing in the coords hash to the output stream
#Is helper function for 
#	The label is modified with respect to the label(id) flag
sub copySubfamSequences {
	(my $dbFile, my $faOut, my $coords_ref, my $idFlag, my $c2sf) = @_;
	my $inseq = Bio::SeqIO->new( -file => $dbFile,    -format => 'fasta' );
	while ( my $seq = $inseq->next_seq ){
		my $id = $seq->display_id(); 
		(my $coords) = $id =~ /^\d+=[^=]+=(\S+:\d+-\d+[+-])/; #get coords from id
		if (exists $coords_ref->{$coords}){ #write the sequence to output if the sequence's coords exist in the input file
			if($c2sf){ #replace subfam from defline with subfam from coordsToDefline hash
				$id =~ s/[^=]+$//; 
				$id .= $c2sf{$coords}; 
			}
			#modify label, unless flag == 0
			if ($idFlag == 1){ #label of "coords=subfam"
				(my $subfam) = $id =~ /\S*=(\S+)$/;
				$seq->display_id($coords ."=". $subfam); 
			}
			elsif ($idFlag == 2){ #label of "coords"
				$seq->display_id($coords); 
			}
			#write seq to output
			$faOut->write_seq($seq);
		}
	}
}
