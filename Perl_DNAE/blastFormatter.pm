#Update 19/8/15 - changed to read from .gz suffix, which is what should be happening. 

use strict; 
use File::Copy;

package blastFormatter; 

sub delZeroIdentityHSPs{
	my $blast = shift; 
	system("gunzip $blast"); 
	my $blastNoZip = $blast; 
	$blastNoZip =~ s/\.gz$//;
	# open(BLAST, "gunzip -c $blast |") or die "open gzipped $blast\n"; #reading with pipe doesn't work with "seek" command...
	open (BLAST, $blastNoZip) || die "couldn't open $blastNoZip\n"; 
	(my $suffix) = $blast =~ /(Seq_\S+)\.gz/; 
	my $tempFile = "temp_$$".$suffix.".txt"; 
	open (TEMP, ">$tempFile") || die "couldn't open $tempFile\n"; 

	my $identityLine; 
	my $identities;  

	while (my $line = <BLAST>){
		if ($line =~ /^\s*Score/){ #new HSP
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
	
	system("gzip -c $tempFile > $blast"); 
	unlink $tempFile; 
	unlink $blastNoZip; 
}

1; 