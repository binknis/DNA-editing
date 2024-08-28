# DNA editing detection in retroelement sequences of a reference genome

If you use the code in this repository, please cite:
Knisbacher BA, Levanon EY. DNA editing of LTR retrotransposons reveals the impact of APOBECs on vertebrate genomes. Molecular Biology and Evolution 2016, 33(2), 554–67.  https://doi.org/10.1093/molbev/msv239

The core pipeline was adapted from code initially provided by Dr. Shai Carmi. 
Thus, please cite also: Carmi S, Church GM, and Levanon EY. 2011. “Large-Scale DNA Editing of Retrotransposons Accelerates Mammalian Genome Evolution.” Nature Communications 2 (November):519.

## Listed below are particularly important scripts used for DNA editing analysis in reference genomes.
Find the files by their name (paths omitted)

#### Genome + UCSC RMSK data preprocessing script
sortGenome.pl

#### Core pipeline scripts
FindClustersByLength.pm

ProcessAlignmentByLength.pm

AnalyzeBlastByLength.pm

runClusterFinder.pl

#### Important helper code
blastFormatter.pm

getAll.pm

#### Important downstream analysis scripts
analysisSubs.pm

createTrackFiles2.pl

createClusterStats3.pl
