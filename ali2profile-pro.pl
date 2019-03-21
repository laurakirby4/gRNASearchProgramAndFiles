#!/usr/bin/perl -w
#use strict;

# read mRNA fasta file and alignment file
# for each fully edited mRNA sequence, count how many alignments cover each position
# the score and length distribution
# only allow one edited mRNA file in the alignment

# 12/14: add the read number; output the best alignment for each position
# 12/16: the fasta file name has a changed format; add the number of reads with the highest score

(@ARGV == 4)|| die "Usage: <alignment file> <mRNA file> <gRNA file> <best alignment file>\n";
my ($alifile, $mRNAfile, $gRNAfile, $bestalifile)=@ARGV;

my $name;
my $fasta="";
my @profile;
my @maxS;
my @minS;
my @maxL;
my @minL;
my @readN;
my @maxSindex;

my @maxreadN;

open(IN, $mRNAfile)||die "cannot open $mRNAfile";
while(<IN>)
{
    chomp;
    
    if($_=~/^>/)
    {
	my @fields=split(/\s+/,$_);
	$name=substr($fields[0],1);
    }
    else
    {
	$fasta.=$_;
    }
}
close(IN);

my %gRNA;
my $gname="";
open(IN, $gRNAfile)||die "cannot open $gRNAfile";
while(<IN>)
{
    chomp;
    if($_=~/^>/)
    {
	my @fields=split(/\s+/,$_);
	$gname=substr($fields[0],1);
    }
    else
    {
	$gRNA{$gname}.=$_;
    }
}
close(IN);

my $size = length($fasta);
for(my $i=0; $i<$size; $i++)
{
    $profile[$i]=0;
    $minS[$i]=999;
    $maxSindex[$i]=-1;
    $maxS[$i]=0;
    $minL[$i]=999;
    $maxL[$i]=0;
    $readN[$i]=0;
    $maxreadN[$i]=0;
}

my @aliarray;
my $index=0;
open(IN, $alifile)||die "cannot open $alifile";
while(my $line=<IN>)
{
    my $mRNA;
    my $record="";
    if($line=~/^>/)
    {
	$record.=$line;
	my $totalread=0;
	my @items = split(/\s+/,$line);
	my @sub_items = split(/:/,$items[2]);
	my $alireads=$sub_items[0];
	$totalread+=$alireads;

	#for(my $k=0; $k<=$#sub_items; $k++)
	#{
	#    my @pairs=split(/\-/,$sub_items[$k]);
	    #print join(" ",@pairs),"\n";
	    #($#pairs==1)||die "$line\n";
	#    if($#pairs==1){
	#	$totalread+=$pairs[1];}
	#    else
	#    {
	#	print STDERR $line;
	#    }
	#}

	#read the next line
	my $nextline=<IN>;
	$record.=$nextline;
	my @fields = split(/\s+/,$nextline);
	#print $fields[7]," x ",$fields[9],"\n";

	for(my $i=$fields[1]; $i<=$fields[2]; $i++)
	{
	    $profile[$i]++;
	    $readN[$i]+=$totalread;

	    if($maxS[$i]<$fields[7])
	    {
		$maxS[$i]=$fields[7];
		$maxSindex[$i]=$index;
		$maxreadN[$i]=$alireads;
		#print "maxS[",$i,"] is ", $maxS[$i],"\n";
	    }
	    elsif($maxS[$i] == $fields[7])
	    {
		$maxreadN[$i]+=$alireads;
	    }
	    
	    if($maxL[$i]<$fields[9])
	    {$maxL[$i]=$fields[9];}

	    if($minS[$i]>$fields[7])
	    {$minS[$i] = $fields[7];}

	    if($minL[$i]>$fields[9])
	    {$minL[$i] = $fields[9];}
	}
	$index++;
	$record.=<IN>;
	$record.=<IN>;
	$record.=<IN>;
    }
    push(@aliarray, $record);
}
close(IN);

($index == ($#aliarray+1))||die "number of alignments is not $index";

open(OUT, ">$bestalifile")||die "cannot create $bestalifile";
#print the best alignment
if($maxSindex[0]!=-1)
{
    print OUT $aliarray[$maxSindex[0]]; 
    my @namelist=split(/\s+/,$aliarray[$maxSindex[0]]);
    print OUT $gRNA{$namelist[2]}, "\n";
}
for(my $i=1; $i<=$#maxSindex; $i++)
{
    if($maxSindex[$i]==$maxSindex[$i-1] or 
       $maxSindex[$i]==-1)
    {
    }
    else
    {
	print OUT $aliarray[$maxSindex[$i]];
	my @namelist=split(/\s+/,$aliarray[$maxSindex[$i]]);
	print OUT $gRNA{$namelist[2]},"\n";
    }
}
close(OUT);

print "> ",$name,"\n";
my @charlist=split(//,$fasta);
for(my $j=0; $j<=$#charlist; $j++)
{
    print $j," ";
    print $charlist[$j]," ";
    print $readN[$j]," ";
    print $maxreadN[$j]," ";
    #print $profile[$j]," ";
    print $maxS[$j]," ";
    print $maxL[$j];
    #print $maxSindex[$j]," ";
    #print $minS[$j]," ";
   
    #print $minL[$j]," ";
    print "\n"
}
print "//";
