#!/usr/bin/env perl

use strict;
use warnings;

if ((!$ARGV[0]) or (!$ARGV[1]) or (!$ARGV[2])){
	print STDERR "\nThis script finds gene IDs (provided in a tab delimited file 1 in column 2) in a tab delimited file (file 2) linking gene ids to some properties and outputs to STDOUT a tab delimited list linking column 1 of file 1 with the user specified column of file 2.

                Usage: prepare_tree_stats_vs_gene_stats.pl <file 1 [gene ID in column2]> <introns_per_gene_summary.txt [gene ID in column 1]> <number of column to be included>\n\n";
        exit 1;
}


my $list_fh=&read_fh ($ARGV[0]);
my $summary_fh=&read_fh ($ARGV[1]);
my $column_number = $ARGV[2];
my (@list_array,@summary_array);
my $original;
my @gene_id;

while(<$list_fh>){
	chomp;
	push(@list_array,$_);
}
while(<$summary_fh>){
	chomp;
	push(@summary_array,$_);
}
#open (OUT,">tree_metrics_vs_makerID.tab") or die $!;



#for(@summary_array){
for (@list_array){
#	$original=$_;
#	my @number=split("\t");
#	my $id=$number[-1];
#	my @id=split("_",$number[-1]);	
#	for (@list_array){
	my @list = split ("\t");
	for (@summary_array){
		if ($_ =~ /^$list[1]/){
			@gene_id=split("\t");
			print "$list[0]\t$gene_id[$column_number-1]\n";
		}
	}
}
#close OUT;

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
sub sum {
	my @array=@_;
	my $sum=0;
	for (@array){
		$sum += $_;
	}
	return $sum;
}
