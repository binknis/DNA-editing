use strict; 

(my $org, my $class, my $fam, my $name) = @ARGV; 

my $seq_count = <STDIN>;
chomp $seq_count; 

my $file = "Name_seq_count_all_orgs.txt"; 
open (OUT, ">>$file") || die "couldn't open $file\n"; 
# print OUT "$org\t$class\t$fam\t$seq_count\n"; 
print OUT "$org\t$class\t$fam\t$name\t$seq_count\n"; 
close(OUT); 