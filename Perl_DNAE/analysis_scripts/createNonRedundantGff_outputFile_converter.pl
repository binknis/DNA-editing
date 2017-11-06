#Function: read the output of a combined list of redundant and non-redundant edited-site and edited-sequence counts and combine to tabular format
#The original file was created by multiple outputs from createNonRedundantGff.pl
#originally written for Nika L2 project (not necessarilly needed elsewhere).
use strict;
use lib "$ENV{HOME}/Perl_DNAE"; 
use getAll;
(my $redundancy_file, my $outfile_name) = @ARGV;
my %seqs = (); 
my %sites = (); 
my %seqs_r = (); 
my %sites_r = ();
my %entries = (); 

my $mut_norm = "GtoA"; 
my $mut_ctrl = "CtoT"; 

#Read input 
my $lines = getAll::lines($redundancy_file); 
die "no lines retrieved for $redundancy_file\n" unless ($lines); 
foreach my $line(@$lines){
	chomp $line; 
	next if $line =~ /^\#/; 
	$line =~ s/1\.00E-/1e-/; #convert excel 1.00E-8 format to 1e-8 format. 
	(my $red, my $mut, my $org, my $class, my $fam, my $pval, my $th, my $site_count, my $seq_count) = split ('\t', $line);
	my $taxa = "$org\t$class\t$fam"; 
	my $params = "$pval\t$th"; 
	if ($red eq 'Non-redundant'){
		$seqs{$mut}{$taxa}{$params} = $seq_count; 
		$sites{$mut}{$taxa}{$params} = $site_count; 
	}
	else {
		$seqs_r{$mut}{$taxa}{$params} = $seq_count; 
		$sites_r{$mut}{$taxa}{$params} = $site_count; 
	}
	$entries{$taxa}{$params} = 0; 
}

#Print output
#Fields: (1)Organism	(2)Class	(3)Family	(4)Pval	(5)Threshold	
#			(6)Edited_seqs_nr	(7)Edited_seqs_nr_control	(8)FP_edited_seqs	(9)Edited_sites_nr	(10)Edited_sites_nr_control	(11)FP_edited_sites
open (OUT, ">$outfile_name") || die "open $outfile_name\n"; 
foreach my $taxa(sort keys %entries){
	foreach my $params(sort keys %{$entries{$taxa}}){
		my $seqs_count = (exists $seqs{$mut_norm}{$taxa}{$params}) ? $seqs{$mut_norm}{$taxa}{$params} : 0; 
		my $sites_count = (exists $sites{$mut_norm}{$taxa}{$params}) ? $sites{$mut_norm}{$taxa}{$params} : 0; 
		my $seqs_ctrl_count = (exists $seqs{$mut_ctrl}{$taxa}{$params}) ? $seqs{$mut_ctrl}{$taxa}{$params} : 0; 
		my $sites_ctrl_count = (exists $sites{$mut_ctrl}{$taxa}{$params}) ? $sites{$mut_ctrl}{$taxa}{$params} : 0; 
		
		my $FP_seqs = ($seqs_count != 0 ) ? ($seqs_ctrl_count / $seqs_count ) * 100 : 0; 
		my $FP_sites = ($sites_count != 0 ) ? ($sites_ctrl_count / $sites_count ) * 100 : 0;
		
		print OUT $taxa."\t".$params."\t";
		print OUT $seqs_count."\t".$seqs_ctrl_count."\t"; 
		printf (OUT "%.3f", "$FP_seqs"); 
		print OUT "\t".$sites_count."\t".$sites_ctrl_count."\t";
		printf (OUT "%.3f", "$FP_sites");
		print OUT "\n"; 			
	}
}
close(OUT); 