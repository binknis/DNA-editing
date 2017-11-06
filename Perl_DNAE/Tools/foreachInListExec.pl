#Function: Execute cmd for each tuple in a file
#			The tuple contains: Organism, Class, Family (mandatory);
#								AND Name, Pvalue (with or w/o '1e-'), Threshold (optional).
#cmd: contains templates that will be replaced with info of each tuple (any subset of templates is OK). 
# ooo - organism; ccc - class; fff - family; nnn - name; ppp - pvalue; ttt - threshold;
#multiple use of any variable or order changes are fine.
use strict; 
use lib "$ENV{HOME}/Perl_DNAE"; 
use getAll; 
(my $listFile, my $cmd) = @ARGV; 

my $list = getAll::lines($listFile); 
die "get list lines failed for $listFile\n" unless $list; 

foreach my $line (@$list){
	my $org, my $class, my $fam, my $name, my $pval, my $th; 
	chomp $line; 
	#parse line data
	#check if lines contain a name field
	my $containsName = 1; 
	my @fields = split(/\s+/,$line);
	if ($#fields < 3 || $fields[3] =~ /^(1e-)?\d+$/){ #doesn't contain name field (list only contains 3 fields or 4th field is the pval).
		$containsName = 0; 
	}
	if ($containsName){
		($org, $class, $fam, $name, $pval, $th) = split(/\s+/,$line); 
	}
	else{
		($org, $class, $fam, $pval, $th) = split(/\s+/,$line);
	}
	$pval = "1e-".$pval unless $pval =~ /1e-/; #for compatability with or w/o '1e-' prefix in pval
	#Replace all tuples requested to be replaced (as specified by templates in cmd)
	my $new_cmd = $cmd; 
	$new_cmd =~ s/ooo/$org/g;
	$new_cmd =~ s/ccc/$class/g;
	$new_cmd =~ s/fff/$fam/g;
	$new_cmd =~ s/nnn/$name/g;
	$new_cmd =~ s/ppp/$pval/g;
	$new_cmd =~ s/ttt/$th/g;
	system($new_cmd);
}
