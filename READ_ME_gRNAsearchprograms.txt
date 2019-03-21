READ ME

fasta-align-pro

This program was written by Dr. Yanni Sun.  This program searches gRNA databases for gRNAs complementary to a given mRNA sequence.  It does not allow mismatches or gaps, but does tolerated G:U base pairs.  This program requires two input parameters to determine search stringency: minimum score and minimum alignment length/gRNA length ratio.  

The score of an alignment is based on the longest stretch of consecutive base pairs with Watson-Crick base pairs scored as 2 points and G:U base pairs scored as 1 point.  For a typical gRNA search, we recommend a minimum score of 45.  For a reduced stringency search we recommend a minimum score of 30.

The minimum alignment length/gRNA length ratio measures the length of the gRNA and mRNA alignment and divides it by the length of the gRNA (excluding the poly-U tail).  For a typical gRNA search, we recommend 0.7 as the minimum ration, but for a reduced stringency search we recommend 0.5.  

To use this program first compile the program using a software compiler such as g++
Example command:
g++ -o fasta-align-pro fasta-align-pro.cc fasta.c nw.c

To perform a gRNA search use the following command structure:

fasta-align-pro gRNAdatabasefile.fa inputmRNAfile.fa minimumalignmentscore minimumlengthratio outputfile.ali.gRNA mismatchespermitted > ouputfile.ali

This command generates two output files:  outputfile.ali.gRNA and outputfile.ali
The first file shows all gRNAs that met the minimum search requirements and the seconds shows the alignments of those gRNAs to the mRNAs

To convert these files into a CSV file, use the following commands:

perl -w ali2profile-pro.pl outputfile.ali inputmRNAfile.fa gRNAdatabasefile.fa outputfile.bestali > outputfile.profile
perl ./data_extraction-v4.pl outputfile.ali outputfile.ali.gRNA outputfile.bestali

This will yield several new files.  outputfile.bestali shows the gRNAs with the best alignments to the mRNA.  outputfile.profile shows the number of gRNA reads covering each nucleotide of the input mRNA sequence.

The csv file ending in (1).csv contains all gRNAs that met the minimum search requirements.  In column A is the first position of the alignment respective to the mRNA, and column B is the last position of the alignment.  Column C is the 5' nt of the gRNA that did not align to the mRNA, column D is the segment of the gRNA that did align to the mRNA, column E is 3' nt of the gRNA that did not align to the mRNA.  Column F is the number of reads for that sequence.  Column G is the unique accession number of that sequence.

The other csv files are groupings based on editing positions and do not contain all edited gRNAs.  

ali2profile-pro.pl was written by Dr. Yanni Sun.
data_extraction-v4.pl was written by Scooter Nowak.  
