#!/usr/bin/env perl

use strict;
use warnings;

if ((!$ARGV[0]) or (!$ARGV[1]) or (!$ARGV[2]) or (!$ARGV[3])){
	print STDERR "\nThis script intersects tree metrics contained in two file that share a gene id in column 2.

                Usage: intersect_tree_metrics.pl <tree_metrics1_vs_gene_id> <tree_metrics2_vs_gene_id> <metric1> <metric2>\n\n";
        exit 1;
}


my $metrics_1_fh=&read_fh ($ARGV[0]);
my $metrics_2_fh=&read_fh ($ARGV[1]);
my $metrics_1=$ARGV[2];
my $metrics_2=$ARGV[3];
#my $taxon_id="Gsa";
my (@metrics_1_array,@metrics_2_array,@line);
my $original;
my @gene_id;

while(<$metrics_1_fh>){
	chomp;
	push(@metrics_1_array,$_);
}
while(<$metrics_2_fh>){
	chomp;
	push(@metrics_2_array,$_);
}
open (OUT_1,">$metrics_1") or die $!;
open (OUT_2,">$metrics_2") or die $!;
open (OUT_csv,">$metrics_1\_vs_$metrics_2.csv") or die $!;

for(@metrics_1_array){
	@gene_id=split("\t");
#	$original=$_;
#	my @number=split("\t");
#	my $id=$number[-1];
#	my @id=split("_",$number[-1]);	
	for (@metrics_2_array){
		if ($_ =~ /$gene_id[1]/){
			@line=split("\t");
			print OUT_1 $gene_id[0]."\n";
			print OUT_2 $line[0]."\n";
			print OUT_csv $gene_id[0] . ",".$line[0]."\n";
		}
	}
}
close OUT_1;
close OUT_2;
close OUT_csv;

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
