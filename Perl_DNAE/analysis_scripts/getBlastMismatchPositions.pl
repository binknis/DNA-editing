#Function: selects a specific list of BLAST pairs froma blast output file and writes the subject-positions of mismatches found in the alignments
#	Typical use - input for R file that plots mismatches (editing vs. other) in alignments. 
#	Notes: 	1. input seqfiles can contain seqs from multiple assemblies
#			2. seqfile fasta IDs must obey original ID from "Data" dir, '='-delimited

use strict; 
use getAll; 
use lib $ENV{HOME} ."/Perl_DNAE/analysis_scripts"; 
use analysisSubs;
use Bio::SeqIO;
use Bio::SearchIO;
use Data::Dumper;

(my $pairsFile, my $blastFile) = @ARGV; 
my @nucs = qw/a c g t/; 
my $printGaps = 0; 

#Get pairs of elements to save their 
my %pairs = (); 
open(PAIRS, $pairsFile) || die "open $pairsFile\n"; 
while (my $l = <PAIRS>){
	chomp $l; 
	(my $assembly, my $source, my $edited) = split(/\t/, $l);
	$pairs{$assembly.",".$source}{$edited}=1; 
}
close(PAIRS);
# print Dumper(\%pairs);
# exit; 
print join("\t", "assembly", "query", "subject", "mm", "pos", "pairInd", "hspCount") ."\n"; 
my $pairIndex = 0; 
my $in = new Bio::SearchIO( -format => 'blast', -file => $blastFile );
while ( my $result = $in->next_result ){
	my $qname = $result->query_name();
	(my $qassembly, my $qcoords) = $qname =~ /^\d+=([^=]+)=([^=]+:\d+\-\d+[+-])=/; #pairs-hash has only coords, and only they will be printed to output
	while ( my $hit = $result->next_hit ){
		my $hitname = $hit->name; 
		(my $hitassembly, my $hitcoords) = $hitname =~ /^\d+=([^=]+)=([^=]+:\d+\-\d+[+-])=/;
		next unless ($qassembly eq $hitassembly and exists $pairs{$qassembly.",".$qcoords}{$hitcoords}); #pair is an editing pair
		$pairIndex++; 
		my $hspCount = 0; 
		while ( my $hsp = $hit->next_hsp ){
			$hspCount++; 
			(my $mm_q, my $mm_s, my $ind, my $q2gaps, my $s2gaps) = analysisSubs::alignment2mmArray(\$hsp);
			foreach my $from(@nucs){
				foreach my $to(@nucs){
					if(exists $mm_s->{$from}{$to}){
						my $mm = uc($to.$from); #NOTE: order is reversed because the mm_array considers the sequence it belongs to as reference.
						foreach my $pos (@{$mm_s->{$from}{$to}}){
							print join("\t", $qassembly, $qcoords, $hitcoords, $mm, $pos, $pairIndex, $hspCount) . "\n";
						}
					}
				}
			}
			if($printGaps){ #print gaps
				foreach my $nuc(@nucs){
					next unless (exists $s2gaps->{$nuc}); 
					my $mm = uc($nuc)."-"; 
					foreach my $pos (@{$s2gaps->{$nuc}}){
						print join("\t", $qassembly, $qcoords, $hitcoords, $mm, $pos, $pairIndex, $hspCount) . "\n";
					}
				}
			}
		}
	}
}
