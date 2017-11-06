use strict; 
use File::Copy;


(my $blast) = @ARGV; 
open (BLAST, $blast) || die "couldn't open $blast\n"; 
(my $suffix) = $blast =~ /(Seq_\S+)/; 
my $tempFile = "temp_$$".$suffix.".txt"; 
open (TEMP, ">$tempFile") || die "couldn't open $tempFile\n"; 

my $identityLine; 
my $identities;  

while (my $line = <BLAST>){
	if ($line =~ /^\s+Score/){ #new HSP
		$identityLine = <BLAST>; 
		$identityLine =~ /Identities = (\d+)/; 
		$identities = $1;
		if ($identities == 0){ #skip this HSP
			while($line = <BLAST>){
				if ($line =~ /^\s*(BLASTN|>|Score)/){
					seek(BLAST,-length($line),1);
					last; 
				}
			}
		}
		else{
			print TEMP $line . $identityLine; 
		}
	}
	else{ #not new HSP
		print TEMP $line; 
	}
}
close(TEMP);
close(BLAST);
File::Copy::move($tempFile, $blast) || die "Didn't move $tempFile back to $blast\n"; 