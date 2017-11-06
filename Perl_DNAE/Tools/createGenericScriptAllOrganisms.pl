use strict; 
use getAll; 
#Notes: filename inputed is without ".sh"
#outname can be a full path

(my $outName, my $numPerFile, my $original_command) = @ARGV; 
#read all organisms from input (whitespace seperated, can be in multiple lines). 

my $organisms = getAll::organisms();  

die "command doesn't have xxx template\n" if $original_command !~ /xxx/; 
my $num = 0; 
my $counter = 0;
 for my $org (@$organisms){
	my $command = $original_command; 
	$command =~ s/xxx/$org/g;
	$counter++; 
	$num++ if ($counter % $numPerFile == 1);
	open (SCRIPT, ">>$outName$num.sh") || die "$outName$num.sh didn't open\n";
	print SCRIPT $command ."\n";
	close (SCRIPT);
}

for (my $i=1; $i<=$num; $i++){
	chmod 0775, "$outName$i.sh"; 
}
