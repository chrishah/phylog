#!/usr/bin/env perl

use strict;
use warnings;

if ((!$ARGV[0]) or (!$ARGV[1])){
	print STDERR "\nThis script intersects tree metrics with intron/exon GC content and intron/exon length.

                Usage: create_tree_metrics_vs_xon_GC_and_length.pl <tree_metrics_vs_gene_id> <gene_id_vs_xon_GC_and_length>\n\n";
        exit 1;
}


my $tree_vs_geneID_fh=&read_fh ($ARGV[0]);
my $gene_vs_intron_GC_length_fh=&read_fh ($ARGV[1]);
#my $taxon_id="Gsa";
my (@tree_vs_geneID_array,@gene_vs_intron_GC_length_array,@line);
my $original;
my @gene_id;

while(<$tree_vs_geneID_fh>){
	chomp;
	push(@tree_vs_geneID_array,$_);
}
while(<$gene_vs_intron_GC_length_fh>){
	chomp;
	push(@gene_vs_intron_GC_length_array,$_);
}
open (OUT_GC,">tree_metric_vs_xon_GC.csv") or die $!;
open (OUT_length,">tree_metric_vs_xon_length.csv") or die $!;

for(@tree_vs_geneID_array){
	@gene_id=split("\t");
#	$original=$_;
#	my @number=split("\t");
#	my $id=$number[-1];
#	my @id=split("_",$number[-1]);	
	for (@gene_vs_intron_GC_length_array){
		if ($_ =~ /$gene_id[1]/){
			@line=split("\t");
			print OUT_GC $gene_id[0].",".$line[1]."\n";
			print OUT_length $gene_id[0].",".$line[2]."\n";
		}
	}
}
close OUT_GC;
close OUT_length;

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
