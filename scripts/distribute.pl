#!/usr/bin/perl
#
#This script distributes files (pattern defined in the script by $selection) into subdirectories of size defined by the user
#

use strict;
use warnings;
use File::Copy qw(move);
use POSIX qw(strftime);
use Cwd qw(abs_path);

my $selection = "fasta"; #alternative would be "aln\.fasta"


my $fasta_dir = $ARGV[0];
my $number_files_per_batch = $ARGV[1];
my $USAGE = "\nUSAGE: ./distribute.pl <PATH/to/DIR> <#files-per-batch>\n";
my @FILES;
my $filecount = 10000;
my @filelist;
 
print strftime("%b %e %H:%M:%S", localtime) . "\n\n";

if (!$ARGV[0]){
        print "$USAGE\n";
        exit;
}else{
	$fasta_dir=abs_path($ARGV[0]);
	print "looking for files ending in \"$selection\" at directory:\n$fasta_dir\n"; 
}

if (!$ARGV[1]){
	$number_files_per_batch = 100;
	print "number of files per batch set to default: $number_files_per_batch\n\n";
}else{
	print "number of files per batch set by user to: $number_files_per_batch\n\n"; 
}

opendir (DIR, "$fasta_dir") or die $!;
@FILES = readdir DIR;
closedir DIR;
chdir "$fasta_dir" or die $!;
for (@FILES){
	chomp;
	if ($_ =~ /$selection$/){
		push (@filelist, $_);
#		print scalar @filelist . "\n";
#		my $remainder = scalar @filelist % $number_files_per_batch;
#		print "remainder $remainder\n";
		if ((scalar @filelist % $number_files_per_batch) == 0){
#			print "create file\n";
			mkdir "$filecount" or die $!;
			chdir "$filecount" or die $!;
			open (FILES,">$filecount.list");
			for (@filelist){
#				print $_ . "\n";
				print FILES $_ . "\n";
				move ("../$_", ".") or die $!;
			}
			close FILES;
			undef @filelist;
			$filecount++;
			chdir "../" or die $!;
		}
	}	
}
if (scalar @filelist >= 1){
	mkdir "$filecount" or die $!;
	chdir "$filecount" or die $!;
	open (FILES,">$filecount");
	for (@filelist){
#		print $_ . "\n";
		print FILES $_ . "\n";
		move ("../$_", ".") or die $!;
	}
	close FILES;
	chdir "../" or die $!;
}
