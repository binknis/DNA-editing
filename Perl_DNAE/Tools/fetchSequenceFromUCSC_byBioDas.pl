use strict;
use Bio::Das;

#### USAGE
unless(@ARGV == 2) {
   die("USAGE: $0 | Input BED file | Output FASTA file\n\n");
}

#### Access files
open(INPUT, "<$ARGV[0]") or die("Could not open input!\n");
open(OUTPUT, ">$ARGV[1]") or die("Could not open output!\n");

#### DAS db
my $das = Bio::Das->new(-source => 'http://genome.cse.ucsc.edu/cgi-bin/das/', -dsn=>'hg19');

#### Line counter
my $count = 0;

#### Work thru each line of BED file
while(defined(my $line = <INPUT>)) {
   #### Make sure line begins with chr
   unless($line =~ /^chr/) { next }

   #### Count lines
   $count++;
   print "#$count $line\n";

   #### Split line
   my ($chr, $start, $end, $label) = '';
   my @line_bits = ();

   @line_bits = split(/\t/, $line);

   $chr = $line_bits[0];
   $start = $line_bits[1];
   $end = $line_bits[2];
   $label = $line_bits[3];
   chomp($label);

   #### Define segment of genome
   my $segment = $das->segment("$chr\:$start\,$end");

   #### Get DNA
   my $dna = $segment->dna;

   #### Print sequence to output
   print OUTPUT "\>$chr:$start-$end\_$label\n$dna\n";
}

#### Close files
close(INPUT);
close(OUTPUT);

exit;