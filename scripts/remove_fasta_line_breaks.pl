#!/usr/bin/env perl

use strict;
use warnings;
use Cwd qw(abs_path);
use POSIX qw(ceil);

if ((!$ARGV[0])) {
	print STDERR "\nThis script removes line breaks from fasta files. if the sequence length exceeds a user defined optional threshold, the sequence will be split into sub sequences of this length. length threshold = 0 does nothing. 

                Usage: remove_fasta_line_breaks.pl *.fasta <length threshold> \n\n";
        exit 1;
}

my $fasta_in=abs_path($ARGV[0]);
#print "path to fasta: $fasta_in\n";
&check_ref_length($ARGV[0],$ARGV[1]);


################################

sub read_fh {
	my $filename = shift @_;
	my $filehandle;
	if ($filename =~ /gz$/) {
		open $filehandle, "gunzip -dc $filename |" or die $!;
	}
	else {
		open $filehandle, "<$filename" or die $!;
	}
	return $filehandle;
}

sub check_ref_length{
        my $ref=$_[0];
##        my $output_filename=$_[1];
        my $critical=$_[1];
        my @header;
        my $header_count=0;
        my (@sequence,@temp_output,@final_output);
        my $full_sequence;
	my $ref_fh=&read_fh($ref);
#        open(REF,"<$ref") or die $!;
        while(<$ref_fh>){
                chomp;
                if ($_ =~ /^>/){
                        push(@header,$_);
#                       print "found header:\n$header[$header_count]\n";
                        $header_count++;
#                       print "header count: $header_count\n";
                        if (@sequence){
                                @temp_output=&finalize_sequence($critical,$header[-2],@sequence);
                                for (@temp_output){
                                        push(@final_output,$_);
                                }
                        }
                        undef @sequence;
                }elsif ($_ =~ /[a-zA-Z]/){
#                       print "found sequence:\n$_\n";
                        push(@sequence,$_);
                }
        }
        @temp_output=&finalize_sequence($critical,$header[-1],@sequence);
        for (@temp_output){
                push(@final_output,$_);
        }
#       print "result:\n";
###        open (OUT,">$output_filename") or die $!;
        for(@final_output){
#               print "$_\n";
###                print OUT "$_\n";
                print "$_\n";
        }
#        close REF;
###        close OUT;

}


sub finalize_sequence{
	my $critical=shift(@_);
	my $header=shift(@_);
	my $full_sequence=join("",@_);
	my $factor;
        my @output;
        if (!$critical){
                $factor=0;
        }else{
                $factor=ceil(length($full_sequence)/$critical);
        }
        if ($factor <= 1){
                push(@output,$header);
                push(@output,$full_sequence);
        }else{ #too long
                print "\nsequence length exceeds limit -> will be split into sub-sequences\n";
                $header=substr $header, 1;
                for (my $i=0; $i<$factor; $i++){
#                        unless ((length(substr $full_sequence, $i*$critical, $critical+31)-31)<0){
                        unless ((length(substr $full_sequence, $i*$critical, $critical))<0){
                                push(@output,">sub$i\_" .$header);
#                                push(@output,substr $full_sequence, $i*$critical, $critical+31);
                                push(@output,substr $full_sequence, $i*$critical, $critical);
                        }
                }
        }
        return @output;
}

