#Function: Execute cmd for each organism and class in Data dir. 
#cmd: if "ooo" template exists - replace with each organism
#	  if "ccc" template exists - replace with each class
#	  if "fff" template exists - replace with each family
#	  if "sss" template exists - replace with each subfamily
# $flags: a regex which includes constraints in the format of 'org=Cat|Dog,family=L1|L2," (see regexes which parse the $flags) 
#script to check it works: perl Tools/foreachExec.pl "echo 'ooo ccc fff sss'"
#Note: in $flags - don't forget the comma at the end of the command! 
#Note: '(' or ')' signs in cmd will be replaced
#example:  perl Tools/foreachExec_new.pl "echo ooo ccc fff sss" "org=Human|Orangutan,subfam=L1PBb," | less
#Note: change Regexes in code to split by "|" and use 'eq' operator to decide if to skip a org/class/fam etc. Regexes cause problems when one org/class/fam/subfam is a substring of another
use strict; 
use FindBin;
my $perlDir = "$FindBin::Bin/.."; # locate this script
use lib "$FindBin::Bin/..";  # use the parent directory of analysis_scripts directory from where this is run
use lib "$FindBin::Bin"; #because this file is in analysis_scripts, this includes that directory
use getAll; 
(my $dataDir, my $cmd, my $flags) = @ARGV; 


my $org_regex, my $class_regex, my $family_regex, my $subfam_regex; 
if ($flags ne ""){
	($org_regex) = $flags =~ /org=([^,]+),/;
	($class_regex) = $flags =~ /class=([^,]+),/;
	($family_regex) = $flags =~ /family=([^,]+),/;
	($subfam_regex) = $flags =~ /subfam=([^,]+),/;
}

my $orgList = getAll::organisms($dataDir); 
die "get orgList failed\n" unless $orgList; 

foreach my $org (@$orgList){
	#create new command with org instead of xxx template
	if ($org_regex) {
		next unless ($org =~ /$org_regex/); 
	}
	my $cmd_with_org = $cmd; 
	$cmd_with_org =~ s/ooo/$org/g; 
	
	my $classList = getAll::classes($dataDir, $org); 
	die "get classList for $org failed\n" unless $classList; 
	
	foreach my $class (@$classList){
		if ($class_regex) {
			next unless ($class =~ /$class_regex/); 
		}
		my $cmd_with_class = $cmd_with_org; 
		$cmd_with_class =~ s/ccc/$class/g; 
		
		my $famList = getAll::families($dataDir, $org, $class); 
			die "get famList for $org $class failed\n" unless $famList; 
			
			foreach my $family (@$famList){
				if ($family_regex) {
					next unless ($family =~ /$family_regex/); 
				}
				my $cmd_with_fam = $cmd_with_class; 
				$cmd_with_fam =~ s/fff/$family/g; 
				
				###insert here : get names, move system command to here. add last for nnn
				my $subfamList = getAll::names($dataDir, $org, $class, $family); 
				die "get subfamList for $org $class $family failed\n" unless $subfamList; 
				foreach my $subfam (@$subfamList){
					# print "here for $org $class $family $subfam\n"; 
					if ($subfam_regex) {
						next unless ($subfam =~ /$subfam_regex/); 
					}
					my $cmd_with_subfam = $cmd_with_fam; 
					$cmd_with_subfam =~ s/sss/$subfam/g; 
					$cmd_with_subfam =~ s/\)/\\\)/g; $cmd_with_subfam =~ s/\(/\\\(/g;
					system($cmd_with_subfam); 
					last if ($cmd !~ /sss/); #no template for subfamily in cmd - don't run per each subfam
				}
				last if ($cmd !~ /fff/); #no template for family in cmd - don't run per each fam
			}
		last if ($cmd !~ /ccc/); #no template for family in cmd - don't run per each class
	}
	last if ($cmd !~ /ooo/); #no template for family in cmd - don't run per each org
}
