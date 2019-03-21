#!/usr/bin/perl
use v5.14;
use strict;

# This program will take as inputs two or three files in the format .gRNA, .ali, and .bestali.
# It will combine the files .gRNA and .ali based on common ID numbers.  The output file
# will be sorted based on the best fit data in the .bestali file and produce a comma
# delinated text file for direct import into excel.
# Columns outputed include 3'-5' identifiers, gRNA sequence, and (who knows what that
# last number is.

# Additional functionality and modifications written by Scooter Nowak

# Take three files as arguments
#(@ARGV == 3) || die "Usage: <file1.ali> <file2.gRNA> <file3.zbestali>\n";

my ($file1,$file2,$file3)=@ARGV;
print "$file1\n$file2\n$file3\n";
# Extract portions of the first input file name to create output file names.
my @outputfile=split(/\./,$file1);

# Set the column locations in each file containing the ID number
# It is not likely these will change
my $ID1=2;
my $ID2=1;

# Set some columns extracted for final csv file.
my $COLUMN1=8;
my $COLUMN2=9;
my $COLUMN4=2;
my $COLUMN6=4;
my $COLUMN7=1;

# Create global arrays and hashes
my @prime;
my %data1;
my %data2;

if (scalar (@ARGV) == 3)
{
    print "Three Arguments! $file1, $file2, $file3\n";
    process_ali_file();
    process_grna_file();
    process_bestali_file();
    align_files();
    generate_csv_file();
    delete_temp_files();
}
elsif (scalar (@ARGV) == 2)
{
    print "Two Arguments! $file1, $file2\n";
    process_ali_file();
    process_grna_file();
    align_files();
    generate_csv_file();
    delete_temp_files();
}

# Open the original ali file and create a rearranged ali file
sub process_ali_file {
    open(IN1, "< $file1") || die "Can't open file $file1!!\n";
    open(OUT1, "> file1output") || die "Can't write file - file1output!!\n";
    while (my $line1 = <IN1>)
    {
        $line1 =~ tr/\n/ /;             # Replace all carriage returns with single space
        $line1 =~ s/> /\n/g;            # Replace the > sympol with a carriage return
        $line1 =~ s/:/ /;
        #print $line1;                   # Debugging
        print OUT1 $line1;              # Print the the output file for later processing
    }
    close (OUT1);                       # Close first output file
    close (IN1);                        # Close first input file
    
    # Process the modified ali file
    open(INA,"< file1output")||die "Can't open $file1";
    <INA>;
    while(<INA>)
    {
        chomp($_);                          # Chomp first newline character
        my @fields=split(/\s+/,$_);         # Split input line and assign values to array "fields"
        $fields[-1] = reverse($fields[-1]); # Reverse the last field so 5' and 3' aligns
        $fields[-1] =~ s/-//g;
        my $family=$fields[$ID1];           # Select a specific field from fields array
        #print $family;                     # Debugging
        shift (@fields);                    # Shift one element off array and move everything left
        my $left=join(" ",@fields);         # Rejoin array fields that were split above
        $data1{$family}=$left;              # Fill hash with rejoined array using value in fields[$ID1] as identifier
        #print "$data1{$family}\n";         # Debugging
    }
    close(INA);                         # Close input file
}

# Open the original gRNA file and create a rearranged gRNA file
sub process_grna_file {
    open(IN2, "< $file2") || die "Can't open file $file2!!\n";
    open(OUT2, "> file2output") || die "Can't write file - file2output!!\n";
    while (my $line2 = <IN2>)
    {
        $line2 =~ tr/\n/ /;             # Replace all carriage returns with single space
        $line2 =~ s/>/\n/g;             # Replace the > symbol with a carriage return
        $line2 =~ s/:/ /;
        #print $line2;                   # Debugging
        print OUT2 $line2;              # Print to the output file for later processing
    }
    close (OUT2);                       # Close second output file
    close (IN2);                        # Close second input file

    # Process the modified gRNA file
    open(INB,"< file2output")||die "Can't open file2output";
    <INB>;
    while(<INB>)
    {
        chomp($_);                      # Chomp first newline character
        my @fields=split(/\s+/,$_);     # Split input line and assign values to array "fields"
        my $family=$fields[$ID2];       # Select a specific field from fields array
        #print "$family\n";                 # Debugging
        shift (@fields);                # Shift one element off array move everything left
        my $left=join(" ",@fields);     # Rejoin array fields that were split above
        $data2{$family}=$left;          # Fill hash with rejoined array using value in fields[$ID2] as identifier
        #print "$data2{$family}\n";     # Debugging
    }
    close(INB);                         # Close input file
}

#Open the original bestali file and create a rearranged bestali file
sub process_bestali_file {
    open(IN3, "< $file3") || die "Can't open file $file3!!\n";
    open(OUT3, "> file3output") || die "Can't write file - file3output!!\n";
    while (my $line3 = <IN3>)
    {
        $line3 =~ tr/\n/ /;             # Replace all carriage returns with single space
        $line3 =~ s/> /\n/g;            # Replace the > sympol with a carriage return
        $line3 =~ s/:/ /;
        #print $line3;                  # Debugging
        print OUT3 $line3;              # Print the the output file for later processing
    }
    close (OUT3);                       # Close second output file
    close (IN3);                        # Close second input file

    # Process the modified bestali file
    open(INC,"< file3output")||die "Can't open $file3";
    <INC>;
    while(<INC>)
    {
        chomp($_);                              # Chomp the first newline character
        my @fields=split(/\s+/,$_);             # Split input line and assign values to array "fields"
        push @prime, @fields[4], @fields[5];    # Push fields[5] and fields[6] into an array
        #print "@fields[5], @fields[6]\n";      # Debugging
    }
    close(INC);                                 # Close input file
}

# Align contents of ali and gRNA files based on ID number
sub align_files {
    open(OUTA, "> exceltemp") || die "Can't write file - exceltemp!!\n";
    foreach my $fname (keys %data2)
    {
        if(exists($data1{$fname}))
        {
            print OUTA $fname," ",$data2{$fname}," | ",$data1{$fname},"\n";   # Output to file
        }
    }
    close(OUTA);                                                              # Close output file
}

# Generate comma deliniated file for import into Excel
sub generate_csv_file {
    open(IN4,"< exceltemp") || die "Can't open exceltemp file";
    open(OUTB, "> $outputfile[0](1).csv") || die "Can't write file - $outputfile[0]_all.csv!!\n";
    #print OUTB "3',5',5' Header, RNA, 3' Tail, #\n";   # Used for debugging
    while(<IN4>)
    {
        my $head;                                       # Create variable to contain the header
        my $tail;                                       # Create variable to contain the tail
        chomp($_);                                      # Chomp off first newline character
        my @fields=split(/\s+/,$_);                     # Split input line containing the header and tail
    
        my $head_tail = $fields[$COLUMN4];              # Assign string to head_tail that contains the header-gRNA-tail
        #print "$head_tail\n";                          # Debugging
        $head_tail =~ s/$fields[-1]/./;                 # Substitute a period for the gRNA section of string
        my @head_tail_split=split(/\./,$head_tail);     # Split the header and tail from each other
        if ($head_tail_split[0] eq "")                  # Check to see if header or tail ar blank, assign space if blank to allow for splitting later
        {
            $head=" ";                                  # If header is blank assign a space
        }
        else
        {
            $head=@head_tail_split[0];                  # If header is not blank assign value to head to be printed later
        }
        if ($head_tail_split[1] eq "")                  
        {
            $tail=" ";                                  # If tail is blank assign a space
        }
        else
        {
            $tail=@head_tail_split[1];                  # If tail is not blank assign value to tail to be printed later
        }
        #print "$head.$tail\n";                         # Used for debugging

        my $rna = $fields[$COLUMN4];                    # Extract gRNA string and assign to rna
        $rna =~ s/$head//;                              # Remove header determined above
        $rna =~ s|(.+)$tail|$1|;                        # Remove tail determined above
    
        # Print selected columns as output
        print OUTB "$fields[$COLUMN1], $fields[$COLUMN2], $head, $rna, $tail, $fields[$COLUMN6], $fields[$COLUMN7]\n";
    }
    close(IN4);                                         # Close input file
    close(OUTB);                                        # Close output file
    
    while(@prime)
    {
        my $threeprime=shift(@prime);                   # Extract three prime value from array prime
        my $fiveprime=shift(@prime);                    # Extract five prime value from array prime
        #print "$threeprime, $fiveprime\n";
        open(OUTC, "> $outputfile[0]($threeprime-$fiveprime).csv") || die "Can't write file - $outputfile[0]($threeprime-$fiveprime).csv!!\n";
        open(IN5,"< $outputfile[0](1).csv")||die "Can't open Outputfile";
        #print "$outputfile[0]\n";
        while(<IN5>)
        {
            my @fields=split(/\s+/,$_);         # Split input line and assign values to array "fields"
            #print "$threeprime, $fiveprime\n";  # Print threeprime and five prime to stdout during run for verification of correct run
            if ($fields[0] >=$threeprime && $fields[1] <=$fiveprime)
            {
                print OUTC $_;                  # Sort the (1) file by the threeprime and fiveprime values
            }
        }
        close(IN5);                             # Close input file
        close(OUTC);                            # Close output file
    }
}

sub delete_temp_files{
    unlink ("file1output");                     # Remove temporary file
    unlink ("file2output");                     # Remove temporary file
    unlink ("file3output");                     # Remove temporary file
    unlink ("exceltemp");                       # Remove temporary file
}