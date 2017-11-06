#Function: creates an interval file from an organisms db directory files
#Input: organism name
#output interval file (stdout)

use strict; 
my $home = $ENV{"HOME"}; 
use lib $home ."/Perl_DNAE"; 
use getAll; 
use Bio::SeqIO;
(my $org) = @ARGV; 

my $classes = getAll::classes($org); 
foreach my $class (@$classes){
	my $fams = getAll::families($org, $class); 
	foreach my $family (@$fams){
		my $subfams = getAll::names($org, $class, $family); 
		foreach my $subfam (@$subfams){
			my $subfam_file = $home ."/Data/". $org ."/". $class ."/db/files_".$family."/Seq_".$subfam;  
			my $in = Bio::SeqIO->new( -file => $subfam_file, -format => 'fasta' );
			while (my $seq = $in->next_seq){
				my $id = $seq->id; 
				$id =~ s/\?//g;
				my @id_parts = split('=',$id);  #1=danRer7=chr1:3147061-3147173+=LINE=L2=CR1-10_DR
				(my $chr, my $start, my $end, my $strand) = $id_parts[2] =~ /(\S+):(\d+)-(\d+)([+-])/; 
				(my $sf, my $c, my $f) = ($id_parts[5], $id_parts[3], $id_parts[4]); 
				print $chr ."\t". $start ."\t". $end ."\t". $strand ."\t". $sf ."\t". $c ."\t". $f ."\n"; 
			}			
		}
	}
}
