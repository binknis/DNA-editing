#Function: creates an ANNOVAR database for a list of specified assemblies and creates BED files for gene files (ensGene and refGene)
#Part 1 --- ANNOVAR DB

#Part 2 --- BED files
#Note: 1. Run from inside Alldbs directory - /home/alu/binknis/Analyses/Annovar/Alldbs (see paths used in createBEDfiles section)
use strict; 

###Consts
my $runAnnovarDBs = 0; 
	my $downdb = 1; 
my $createBEDfiles = 1; 

(my $assemblies) = @ARGV; 
$assemblies =~ s/,$//; #enable trailing comma in input
my @assembly_arr = split (',', $assemblies); 

if($runAnnovarDBs){
	foreach my $assembly (@assembly_arr){
		###  Download gene annotation
		#Refseq gene annotation
		print "******************** comm1_refseq for $assembly ************************\n"; 
		my $comm1_refseq = "perl516 /private/common/Software/ANNOVAR/annovar/annotate_variation.pl -downdb -buildver ".$assembly." gene ".$assembly."db";
		system ($comm1_refseq); 
		#Ensembl gene annotation
		print "******************** comm1_ensgene for $assembly ************************\n"; 
		my $comm1_ensgene = "perl516 /private/common/Software/ANNOVAR/annovar/annotate_variation.pl -downdb -buildver ".$assembly." ensgene ".$assembly."db";
		system ($comm1_ensgene); 
		
		### Download genome sequence - single command for all gene formats
		if($downdb){
			print "******************** comm2 for $assembly ************************\n"; 
			my $comm2 = "perl516 /private/common/Software/ANNOVAR/annovar/annotate_variation.pl --buildver ".$assembly." --downdb seq ".$assembly."db/".$assembly."_seq"; 
			system ($comm2); 
		}
		
		### 
		#Refseq-specific command to build the mRNA fasta files from the refGene file
		my $comm3_refseq; 
		if (-e $assembly."db/".$assembly."_seq/".$assembly.".fa"){
			print "******************** comm3 refseq for $assembly (opt1) ************************\n"; 
			$comm3_refseq = "perl516 /private/common/Software/ANNOVAR/annovar/retrieve_seq_from_fasta.pl ".$assembly."db/".$assembly."_refGene.txt -seqfile ".$assembly."db/".$assembly."_seq/".$assembly.".fa -format refGene -outfile ".$assembly."db/".$assembly."_refGeneMrna.fa"; 
		}
		else{
			print "******************** comm3 refseq for $assembly (opt2) ************************\n"; 
			$comm3_refseq = "perl516 /private/common/Software/ANNOVAR/annovar/retrieve_seq_from_fasta.pl ".$assembly."db/".$assembly."_refGene.txt -seqdir ".$assembly."db/".$assembly."_seq -format refGene -outfile ".$assembly."db/".$assembly."_refGeneMrna.fa"; 
		}
		system ($comm3_refseq); 
		
		#ensGene-specific command to build the mRNA fasta files from the refGene file
		my $comm3_ensgene; 
		if (-e $assembly."db/".$assembly."_seq/".$assembly.".fa"){
			print "******************** comm3 ensgene for $assembly (opt1) ************************\n"; 
			$comm3_ensgene = "perl516 /private/common/Software/ANNOVAR/annovar/retrieve_seq_from_fasta.pl ".$assembly."db/".$assembly."_ensGene.txt -seqfile ".$assembly."db/".$assembly."_seq/".$assembly.".fa -format ensGene -outfile ".$assembly."db/".$assembly."_ensGeneMrna.fa"; 
		}
		else{
			print "******************** comm3 ensgene for $assembly (opt2) ************************\n"; 
			$comm3_ensgene = "perl516 /private/common/Software/ANNOVAR/annovar/retrieve_seq_from_fasta.pl ".$assembly."db/".$assembly."_ensGene.txt -seqdir ".$assembly."db/".$assembly."_seq -format ensGene -outfile ".$assembly."db/".$assembly."_ensGeneMrna.fa"; 
		}
		system ($comm3_ensgene); 
		
	}
}

###Download additional files from UCSC into assembly database
#"cytoBand.txt.gz"  - 
#"all_ests" - Expression
#"targetScanS.txt.gz" - predicted miRNA targets

#Create separate files for: Gene (tx regions); CDS (); 
my $num_fields_with_bin = 16; #CONST
my $num_bp_upstream = 1000; #CONST. number of basepairs to fetch when fetching "upstream" coordinates for promoter regions

if ($createBEDfiles){
	foreach my $assembly (@assembly_arr){
		my $beddir = $assembly."db/BEDs";
		mkdir  $beddir; 
		my @genedbs = ("refGene", "ensGene", "xenoRefGene");
		# my @genedbs = ("xenoRefGene");
		foreach my $genedb (@genedbs){
			my $genedb_file = $assembly."db/".$assembly."_".$genedb.".txt"; 
			next unless (-e $genedb_file); 
			open (my $genedb_fh, $genedb_file) or die "open $genedb_file\n"; 
			
			my $mrna_file = $beddir ."/".$genedb."_mrna.bed"; 
			my $cds_file = $beddir ."/".$genedb."_cds.bed";
			my $exon_file = $beddir ."/".$genedb."_exons.bed"; 
			my $intron_file = $beddir ."/".$genedb."_introns.bed"; 
			my $utr5_file = $beddir ."/".$genedb."_5utr.bed"; 
			my $utr3_file = $beddir ."/".$genedb."_3utr.bed"; 
			my $cdsExon_file = $beddir ."/".$genedb."_cdsExons.bed"; 
			my $upstream_file = $beddir ."/".$genedb."_upstream".$num_bp_upstream."bp.bed"; 
			my $tss_file = $beddir ."/".$genedb."_tss.bed"; 
			
			
			open (my $mrna_fh, ">".$mrna_file) or die "open $mrna_file\n"; 
			open (my $cds_fh, ">".$cds_file) or die "open $cds_file\n"; 
			open (my $exon_fh, ">".$exon_file) or die "open $exon_file\n"; 
			open (my $cdsExon_fh, ">".$cdsExon_file) or die "open $cdsExon_file\n"; 
			open (my $intron_fh, ">".$intron_file) or die "open $intron_file\n"; 
			open (my $utr5_fh, ">".$utr5_file) or die "open $utr5_file\n"; 
			open (my $utr3_fh, ">".$utr3_file) or die "open $utr3_file\n"; 
			open (my $upstream_fh, ">".$upstream_file) or die "open $upstream_file\n"; 
			open (my $tss_fh, ">".$tss_file) or die "open $tss_file\n"; 
			
			##Read mRNA input file
			while(my $l = <$genedb_fh>){
				chomp $l; 
				my @f = split(/\t/, $l); 
				die "incompatible amount of fields in $genedb_file - $#f\n" unless ($#f==$num_fields_with_bin-1 or $#f==$num_fields_with_bin-2); 
				shift @f if  ($#f==$num_fields_with_bin-1); #discard "bin" column if exists
				(my $name, my $chrom, my $strand,
					my $txStart, my $txEnd, my $cdsStart, my $cdsEnd, 
					my $exonCount, my $exonStarts, my $exonEnds, 
					my $score, my $name2, my $cdsStartStat, my $cdsEndStat, my $exonFrames) = @f; 
				my $outname = join('|', $name, $name2); 
				
				##Output
				#full mRNA BED
				print $mrna_fh join("\t", $chrom, $txStart, $txEnd, $outname, ".", $strand), "\n";
				#full CDS BED  ; 5' and 3' UTRs
				if ($cdsStart != $cdsEnd){ #avoid output of 'empty' CDSs (the format writes the same position for start and end of CDS when there is no CDS).
					print $cds_fh join("\t", $chrom, $cdsStart, $cdsEnd, $outname, ".", $strand), "\n"; #CDS BED
					if ($txStart != $cdsStart){
						print $utr5_fh join("\t", $chrom, $txStart, $cdsStart, $outname, ".", $strand), "\n"; #5' UTR BED
					}
					if ($cdsEnd != $txEnd){
						print $utr3_fh join("\t", $chrom, $cdsEnd, $txEnd, $outname, ".", $strand), "\n"; #3' UTR BED
					}
					
				}
				#exon BED
				$exonStarts =~ s/,$//; $exonEnds =~ s/,$//; #remove trailing comma
				my @exonStarts = split(',', $exonStarts); 
				my @exonEnds = split(',', $exonEnds); 
				for my $i (0 .. $#exonStarts){
					print $exon_fh join("\t", $chrom, $exonStarts[$i], $exonEnds[$i], $outname, ".", $strand), "\n"; 
					#coding exons BED (UTR regions are stripped from 1st and last exons)
					if ( ($exonStarts[$i] >= $cdsStart and $exonStarts[$i] < $cdsEnd) or ($exonEnds[$i] > $cdsStart and $exonEnds[$i] <= $cdsEnd) ){ #exon overlaps CDS
						my $exStart = ($exonStarts[$i] > $cdsStart ? $exonStarts[$i] : $cdsStart); #get only overlapping region with CDS - truncate 5' UTR
						my $exEnd = ($exonEnds[$i] < $cdsEnd ? $exonEnds[$i] : $cdsEnd); #get only overlapping region with CDS - truncate 3' UTR
						print $cdsExon_fh join("\t", $chrom, $exStart, $exEnd, $outname, ".", $strand), "\n"; #CDS regions in exons
					}
				}
				#intron BED
				for my $i (1 .. $#exonStarts){
					print $intron_fh join("\t", $chrom, $exonEnds[$i-1], $exonStarts[$i], $outname, ".", $strand), "\n"; 
				}
				
				#upstream BED
				if($txStart>0){
					my $promoterStart = $txStart - $num_bp_upstream; 
					$promoterStart = 0 if $promoterStart < 0; 
					print $upstream_fh join("\t", $chrom, $promoterStart, $txStart, $outname, ".", $strand), "\n"; #'upstream' excludes the tss (txStart is 0-base)
				}
				#TSS BED
				print $tss_fh join("\t", $chrom, $txStart, $txStart+1, $outname, ".", $strand), "\n"; #txStars is 0-based so add 1 to 'end'
			}
			
			close($mrna_fh); 
			close($cds_fh); 
			close($exon_fh); 
			close($cdsExon_fh); 
			close($intron_fh); 
			close($utr5_fh); 
			close($utr3_fh); 
			close($upstream_fh); 
			close($tss_fh); 
			close($genedb_fh);
			
			#Sort output files
			foreach my $file ($mrna_file, $cds_file, $exon_file, $cdsExon_file, $intron_file, $utr5_file, $utr3_file, $upstream_file, $tss_file){
				if (-e $file){
					system("sortBed -i $file > $file"."_temp"); 
					system("mv $file"."_temp $file"); 
				}
			}
		}
	}
}
