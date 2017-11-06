#Function: Execute cmd for each organism and class in Data dir. 
#cmd: if "ooo" template exists - replace with each organism
#	  if "ccc" template exists - replace with each class
#script to check it works: perl Tools/foreachExec.pl "echo 'ooo ccc zzz'"
use strict; 
use getAll; 
(my $cmd) = @ARGV; 

my $orgList = getAll::organisms(); 
die "get orgList failed\n" unless $orgList; 

foreach my $org (@$orgList){
	#create new command with org instead of ooo template
	my $cmd_with_org = $cmd; 
	$cmd_with_org =~ s/ooo/$org/g; 
	
	if ($cmd_with_org =~ /ccc/){
		my $classList = getAll::classes($org); 
		die "get classList for $org failed\n" unless $classList; 
	
		foreach my $class (@$classList){
			my $cmd_with_class = $cmd_with_org; 
			$cmd_with_class =~ s/ccc/$class/g; 
			if ($cmd_with_class =~ /zzz/){
				my $famList = getAll::families($org, $class); 
				die "get famList for $org $class failed\n" unless $famList; 
				
				foreach my $fam (@$famList){
					my $cmd_with_fam = $cmd_with_class; 
					$cmd_with_fam =~ s/zzz/$fam/g; 
					system($cmd_with_fam); 
				}
			}
			else{ #cmd for org and class only
				system($cmd_with_class); 
			}
		}
	}
	else{
		system($cmd_with_org);
	}
}
