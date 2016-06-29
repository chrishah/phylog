#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;

if ((!$ARGV[0])){
        print STDERR "\nThis script calculates the the average bootstrap support of all trees in a given directory inferred by RAxML (*bipartitionsBranchLabels* format i.e. bootstrap values in []).

                Usage: calculate_average_bootstrap_support.pl <directory>\n\n";
        exit 1;
}

my (@split_tree1,@bootstrap,@summary);
my $fasta_dir=$ARGV[0];
my $tree_fh;
my $path=getcwd;

################

#my $treefh= &read_fh ($ARGV[1]);

opendir (DIR, "$fasta_dir") or die $!;
my @FILES = readdir DIR;
closedir DIR;
chdir "$fasta_dir" or die $!;
for (@FILES){
        chomp;
        if ($_ =~ /^RAxML/){
		undef @bootstrap;
		$tree_fh=&read_fh ($_);
		my $gene_name=$_;
		print $gene_name . "\n";
		my @temp=split("_");
#		print "$temp[-1]\n";
		my @id=split('\.', $temp[-1]);
#		print "$id[0]\n";
		while(<$tree_fh>){
			chomp;
#			print "$_\n";
			@split_tree1 = split('\[');
		}
		shift(@split_tree1);
		for (@split_tree1){
			if ($_ =~ /([0-9]+)/){
#				print "$1\n";
				push (@bootstrap,$1);
			}
		}
		undef my @cluster_name;
		push (@cluster_name, $temp[-2], $id[0]);
		my $cluster_id = join("_", @cluster_name);
		push (@summary,&average(@bootstrap)."\t".$gene_name."\t".$id[0]."\t".$cluster_id.".");
	}
}
@summary=sort{ $a <=> $b } @summary;
open(OUT,">$path/SUMMARY_average_bootstrap.txt") or die $!;
for (@summary){
        print OUT "$_\n";
}
close OUT;


#while (<$treefh>){
#	chomp;
#	@split_tree1 = split(":");
#}
#shift(@split_tree1);
#
#for(@split_tree1){
#	print "$_\n";
#}
#print "\n\n";
#
#for (@split_tree1){
#	if ($_ =~/([0-9.]+)/){
#		print "$1\n";
#		push (@branch_lengths,$1);
#	}
#}
#
#print "\nthe cummulative branch length is: ". &sum(@branch_lengths) ."\n";
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
sub average {
	my @array=@_;
	my $sum=0;
	for (@array){
		$sum += $_;
	}
	my $average = $sum / scalar @array;
	return $average;
}
