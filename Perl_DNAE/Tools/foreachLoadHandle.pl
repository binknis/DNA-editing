#Function: Execute cmd for each organism and class in Data dir and load handle. 
#cmd: if "ooo" template exists - replace with each organism
#	  if "ccc" template exists - replace with each class
#	  if "fff" template exists - replace with each family
#	  if "sss" template exists - replace with each subfamily
# $flags: a regex which includes constraints in the format of 'org=Cat|Dog,fam=L1|L2," (see regexes which parse the $flags) 
# $processes: amount of processes to keep busy at all times, running the specified commands (default = 12)
#script to check it works: perl Tools/foreachExec.pl "echo 'ooo ccc fff sss'"
#Note: in $flags - don't forget the comma at the end of the command! 
#Note: '(' or ')' signs in cmd will be replaced in loadHandler.pl
#Note: This script finishes and loadHandler.pl is started with nohup
#example:  perl Tools/foreachExec_new.pl "echo ooo ccc fff sss" "org=Human|Orangutan,subfam=L1PBb," | less
use strict; 
use lib "$ENV{HOME}/Perl_DNAE"; 
use getAll; 
(my $cmd, my $flags) = @ARGV; 

my $org_regex, my $class_regex, my $family_regex, my $subfam_regex, my $processes; 
if ($flags ne ""){
	($org_regex) = $flags =~ /org=([^,]+),/;
	($class_regex) = $flags =~ /class=([^,]+),/;
	($family_regex) = $flags =~ /family=([^,]+),/;
	($subfam_regex) = $flags =~ /subfam=([^,]+),/;	
	($processes) = $flags =~ /processes=([^,]+),/;  
}
unless ($processes){ #set default
	$processes = 12;
}


my $tmpFile = "tmp_in_use_dont_erase_".$$.".txt"; 
open (TMP, ">$tmpFile") or die "couldn't open $tmpFile\n"; 
my $orgList = getAll::organisms(); 
die "get orgList failed\n" unless $orgList; 

foreach my $org (@$orgList){
	#create new command with org instead of xxx template
	if ($org_regex) {
		next unless ($org =~ /$org_regex/); 
	}
	my $cmd_with_org = $cmd; 
	$cmd_with_org =~ s/ooo/$org/g; 
	
	my $classList = getAll::classes($org); 
	die "get classList for $org failed\n" unless $classList; 
	
	foreach my $class (@$classList){
		if ($class_regex) {
			next unless ($class =~ /$class_regex/); 
		}
		my $cmd_with_class = $cmd_with_org; 
		$cmd_with_class =~ s/ccc/$class/g; 
		
		my $famList = getAll::families($org, $class); 
			die "get famList for $org $class failed\n" unless $famList; 
			
			foreach my $family (@$famList){
				if ($family_regex) {
					next unless ($family =~ /$family_regex/); 
				}
				my $cmd_with_fam = $cmd_with_class; 
				$cmd_with_fam =~ s/fff/$family/g; 
				
				###insert here : get names, move system command to here. add last for nnn
				my $subfamList = getAll::names($org, $class, $family); 
				die "get subfamList for $org $class $family failed\n" unless $subfamList; 
				foreach my $subfam (@$subfamList){
					# print "here for $org $class $family $subfam\n"; 
					if ($subfam_regex) {
						next unless ($subfam =~ /$subfam_regex/); 
					}
					my $cmd_with_subfam = $cmd_with_fam; 
					$cmd_with_subfam =~ s/sss/$subfam/g; 
					#$cmd_with_subfam =~ s/\)/\\\)/g; $cmd_with_subfam =~ s/\(/\\\(/g; ###this is done in the LOAD HANDLER
					print TMP $cmd_with_subfam ."\n"; 
					last if ($cmd !~ /sss/); #no template for subfamily in cmd - don't run per each subfam
				}
				last if ($cmd !~ /fff/); #no template for family in cmd - don't run per each fam
			}
		last if ($cmd !~ /ccc/); #no template for family in cmd - don't run per each class
	}
	last if ($cmd !~ /ooo/); #no template for family in cmd - don't run per each org
}
close (TMP); 

#run load handler
system("nohup perl Tools/loadHandler.pl $tmpFile $processes > nohup.out.loadHandle" . $$ . " 2>&1"); 
unlink($tmpFile); 
