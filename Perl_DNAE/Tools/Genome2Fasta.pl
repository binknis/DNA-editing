use strict;
use Bio::SeqIO;
( my $genome_file, my $interval_file, my $fasta_out, my $lc_flag) = @ARGV;
my %chr_objs; 
my $sequence;
my $subseq;
my $chr_obj; 
my $seq_obj; 
my $genoName;
my $fasta_id; 
$genome_file =~ /(\S+)_(\S+)\.(\S+)/;
$genoName = $2;
#print "$genoName";  
open (INTERVAL_FILE ,$interval_file) || die ("couldn't open interval file");
$fasta_out = "fasta.txt" if $fasta_out eq "";

#open (my $file ,">$fasta_out");
my $seqio_obj = Bio::SeqIO->new(-file => ">$fasta_out", -format => 'fasta' );

my $inseq = Bio::SeqIO->new( -file => $genome_file,    -format => 'fasta' );
while ( my $seq = $inseq->next_seq() )
{
        my $id = $seq->display_id();
        #print "$id\n";
		$seq_obj = Bio::Seq->new(-seq => $seq->seq(), -alphabet => 'dna' );
        $chr_objs{$id} = $seq_obj; 
}

# my $count=0; 

while(my $line  = <INTERVAL_FILE>)
{
 next if $line =~ /^#/; 
 chomp $line;
 (my $chr, my $genoStart, my $genoEnd , my $strand) = split("\t",$line);
 
	$chr_obj = $chr_objs{$chr};
	$subseq = $chr_obj->subseq($genoStart+1,$genoEnd);

    if($strand  =~ /-/)
    {
     $seq_obj = Bio::Seq->new(-seq => $subseq, -alphabet => 'dna' );
         
         $seq_obj = $seq_obj->revcom();
         
         $subseq = $seq_obj->seq();
    }
	if ($lc_flag){
		$subseq = lc $subseq; 
	}
    #$fasta_id = ">".$genoName."_".$chr."_".$genoStart."_".$genoEnd."_".$strand;
	$fasta_id = $genoName."_".$chr."_".$genoStart."_".$genoEnd."_".$strand;
    #print $file "$fasta_id \n";    
	$seq_obj = Bio::Seq->new(-seq => $subseq, -display_id =>$fasta_id, -alphabet => "dna" );              

	$seqio_obj->write_seq($seq_obj);

	#print $file "$subseq \n"; 
        # $count++; 
        # last if $count > 30; 
}
