#Function: creates shell-scripts to run completeBlasts.pl for all organisms that aren't finished 
#@ARGV: directory in which organisms are stored (typically "../../Data")
#Note: at this point it creates scripts for LINE LTR SINE completion only. (see regex below)

use strict; 
use getAll; 

my $orgs = getAll::organisms(); 
my %incomplete = (); 
my %complete = (); 

foreach my $org (@$orgs){
	my $classes = getAll::classes($org); 
	foreach my $class (@$classes){
		next unless $class =~ /(LINE|LTR|SINE)/; 
		my $fams = getAll::families($org, $class); 
		foreach my $fam (@$fams){
			my $clusts = getAll::clustersFam($org, $class, $fam); 
			if ($clusts == 0){
				$incomplete{$org}{$class} = () unless exists $incomplete{$org}{$class}; 
				push (@{$incomplete{$org}{$class}}, $fam); 
			}
			else{
				$complete{$org}{$class} = () unless exists $complete{$org}{$class}; 
				push (@{$complete{$org}{$class}}, $fam); 
			}
		}
	}
}


print "INCOMPLETE FAMILIES:\n"; 
foreach my $org(sort keys %incomplete){
	foreach my $class (sort keys %{$incomplete{$org}}){
		print "$org $class: @{$incomplete{$org}{$class}}\n"; 
	}
}
print "********************************************************\n"; 
print "COMPLETE FAMILIES:\n";
foreach my $org(sort keys %complete){
	foreach my $class (sort keys %{$complete{$org}}){
		print "$org $class: @{$complete{$org}{$class}}\n"; 
	}
}

#creates scripts to run FamilyRunClusterFinder.pl for families (11 per file):
my $count =0; 
my $fileNum = 1; 

foreach my $org(sort keys %incomplete){
	foreach my $class (sort keys %{$incomplete{$org}}){
		foreach my $fam (sort @{$incomplete{$org}{$class}}){
			my $script	= "script_FamilyRunClusterFinder".$fileNum.".sh"; 
			open (SCRIPT, ">>$script") || print "couldn't open script file\n"; 
			print SCRIPT "nohup perl family_scripts/FamilyRunClusterFinder.pl $org 1 $fam 3 16 3 16 11 $class &\n";
			close (SCRIPT);
			chmod 744, $script; 
			$count++; 
			$fileNum++ if $count % 11 == 0; 
		}
	}
}
