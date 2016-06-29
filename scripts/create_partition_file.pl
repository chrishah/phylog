#!/usr/bin/perl
#This script takes a modified version of the info file returned after FASconCAT concatenation and the summary file created by the monophyly_masking.pl script and produces a partition file to be used with RAxML
use strict;
use warnings;
use POSIX qw(strftime);
use File::Copy qw(copy move);
use Getopt::Long;

#my $fasta_dir = "/home/christoph/my_orthomcl/Inflation-2.1/cluster_selection_results/";
my $FASconCAT_info = $ARGV[0];
my $summary_file = $ARGV[1];
my (@FASconCAT, @summary_file, @summary_line_split, @partition_file);
my ($help, $line);

my $USAGE = "\nThis script takes a modified version of the info file returned after FASconCAT concatenation and the summary file created by the monophyly_masking.pl script and produces a partition file to be used with RAxML

		USAGE: ./create_partition_file.pl <modified FASconCAT-info file> <summary-output tab file>

		options:\n
		--help		show this infomation\n";

GetOptions (	"help!" => \$help) or die "Incorrect usage!\n$USAGE";




if ((!$ARGV[0]) || (!$ARGV[1]) ||  ($help)){
	print "$USAGE\n";
	exit;
}

open (FASconCAT,"<$ARGV[0]") or die $!;
open (SUMMARY,"<$ARGV[1]") or die $!;

map {push(@FASconCAT, $_)} <FASconCAT>;
map {push(@summary_file, $_)} <SUMMARY>;
for (@FASconCAT){
	chomp;
	my $whole_line = $_;
	my @FASconCAT_line = split / /;
	my $FASconCAT_ID = $FASconCAT_line[0] . ".aln";
	print "$FASconCAT_ID ... ";
	for (@summary_file){
		chomp;
		if ($_ =~ $FASconCAT_ID){
#			print "$_\n";
			undef $line;
			undef @summary_line_split;
			@summary_line_split = split /\t/;
#			print "$summary_line_split[3]\n";
			$line = $summary_line_split[3] . ", " . $whole_line;
			push (@partition_file, $line);
#			print "$line\n";
			if (!$summary_line_split[3]){
				print "no model found\n";
				open (ERROR,">>problematic.log") or die $!;
				print ERROR $FASconCAT_ID . "\n";
				close ERROR;
			}else {
				print "$summary_line_split[3]\n";
			}
		}
	}

#	my $summary_line = grep { $_ =~ $FASconCAT_line[0]} @summary_file;
#	print "$summary_line\n";
#	map {chomp; print "$_\n";} @summary_line;
#	my @summary_line_split = split (/\t/, $summary_line);
#	my @summary_line_split = split (/\t/, @summary_line);
#	print "$summary_line_split[3]\n";
}

open (PARTITION,">part.raxml") or die $!;
map {print PARTITION "$_\n"} @partition_file;

close FASconCAT;
close SUMMARY;
