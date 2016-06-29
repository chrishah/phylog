#!/usr/bin/perl
#This script crops the sequence IDs of in every *.phy file in the given directory to a number of characters defined by the user


use strict;
use warnings;
use POSIX qw(strftime);
use File::Copy qw(copy move);
use Getopt::Long;

my $phylip_dir = $ARGV[0];
my $number_of_characters = $ARGV[1];
my (@FILES, @line, @ID, @content, @FILE);
my ($help, $ID, $line, $count);

my $USAGE = "\nThis script crops the sequence IDs of in every *.phy file in the given directory to a number of characters defined by the user

		USAGE: ./crop_IDs.pl <PATH/to/DIR> <number of characters to retain>\n\n";

GetOptions (	"help!" => \$help) or die "Incorrect usage!\n$USAGE";


if ((!$ARGV[0]) || (!$ARGV[1]) || ($help)){
	print "$USAGE\n";
	exit;
}

opendir (DIR, "$phylip_dir") or die $!;
@FILES = readdir DIR;
closedir DIR;

chdir "$phylip_dir" or die $!;
mkdir "cropped" or die $!;
print strftime("%b %e %H:%M:%S", localtime) . "\n\n";
for (@FILES){
        chomp;
        if ($_ =~ /phy$/){
		print "$_ ... ";
		open (PHYLIPFILE,"<$_") or die $!;
		open (CROPPEDPHYLIPFILE,">cropped/$_") or die $!;
		map {chomp; push (@content, $_)} <PHYLIPFILE>;
		for (my $i = 0; $i < @content; $i++){
			chomp;
			if ($i == 0){
				push (@FILE, $content[$i]);	
			}elsif ($i >= 1){
				@line = split (/\t/, $content[$i]);
				$ID = substr $line[0], 0, $ARGV[1];
				$line = $ID . "\t" . $line[1];
				push (@FILE, $line);	
			}
		}	
		close PHYLIPFILE;
		map {print CROPPEDPHYLIPFILE "$_\n"} @FILE;
		undef @FILE;
		undef @content;
		close CROPPEDPHYLIPFILE;
		print "DONE!\n";
		$count++;
	}
}
print "Cropped $count files - FINISHED!\n";
print strftime("%b %e %H:%M:%S", localtime) . "\n";
