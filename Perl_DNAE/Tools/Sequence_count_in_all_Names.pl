#Function: counts the amount of sequences in each Family and subfamily in Data folder
#	Input: none.  
#	Output: 2 files - 	1. Name_seq_count_all_orgs.txt: Number of sequences in each Name.
#						2. Family_seq_count_all_orgs.txt: Number of sequences in each Family.
#	Note: needs Tools/sequence_count_printer.pl script to work (it generates the first output file). 

use strict; 
use getAll; 

#create file with count of sequences in each name (one file for all organisms, classes and families. line for each name)
my $orgs = getAll::organisms(); 
foreach my $org(@$orgs){
	my $classes = getAll::classes($org);
	next if $classes == 0; 
	foreach my $class (@$classes){
		my $fams = getAll::families($org, $class); 
		next if $fams == 0; 
		foreach my $fam (@$fams){
			my $name_files = getAll::name_files($org, $class, $fam, 1); 
			next if $name_files == 0; 
			foreach my $nf (@$name_files){
				(my $name) = ($nf =~ /Seq_(\S+)$/);
				# print "$org\t$class\t$fam\t$name\n"; 
				my $cmd = "tac $nf | grep -m 1 '>' | sed 's/>//' | sed 's/=.*//' | perl Tools/sequence_count_printer.pl $org $class $fam $name"; 
				$cmd =~ s/\(/\\\(/g; $cmd =~ s/\)/\\\)/g; 
				system ($cmd); 
			}
		}
	}
}
### create file with count of sequences in each family (for all classes and organisms - all 1 file!) ###
#read count from name seq count file
my %famCount = (); 
my $name_lines = getAll::lines("Name_seq_count_all_orgs.txt"); 
my $famCount=0;
foreach my $nl (@$name_lines){
	chomp $nl; 
	(my $org, my $class, my $fam, my $name, my $nameCount) = split (/\s+/,$nl); 
	$famCount{$org}{$class}{$fam} = 0 unless exists $famCount{$org}{$class}{$fam};
	$famCount{$org}{$class}{$fam} += $nameCount; 
}

#create fam seq count file
my $fam_seq_count_all_organisms = "Family_seq_count_all_orgs.txt"; 
open (FAMFILE, ">$fam_seq_count_all_organisms") || die "didn't open $fam_seq_count_all_organisms\n"; 
foreach my $org (sort keys %famCount){
	print $org ."\n"; 
	foreach my $class (sort keys %{$famCount{$org}}){
		foreach my $fam (sort keys %{$famCount{$org}{$class}}){
			print FAMFILE $org."\t".$class."\t".$fam."\t".$famCount{$org}{$class}{$fam}."\n";
		}
	}
}
close(FAMFILE);
