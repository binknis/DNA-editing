use strict;
use warnings;
use LWP::Simple;
use XML::XPath;
use XML::XPath::XMLParser;

#Use DAS of UCSC to fetch specific sequence by its given chromosome position

my $URL_gene ="http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr2:100000,100100";

my $xml = get($URL_gene);

my $xp = XML::XPath->new(xml=>$xml);

my $nodeset = $xp->find('/DASDNA/SEQUENCE/DNA/text()'); # find all sequences
# there should be only one node, anyway:    
foreach my $node ($nodeset->get_nodelist) {
   my $seq = $node->getValue;
   $seq =~ s/\s//g; # remove white spaces
   print $seq, "\n";
}