#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);
use Cwd;

my $path=getcwd;
my $list_file;
my $summary_file;
my $taxon;
my (@list_array, @summary_array, @summary_array_split,@out_array,@taxa);
GetOptions (
#	"blast2goannotfile:s" => \$blast2goannotfile,
	"b:s"    => \$list_file,
	"s:s"	=> \$summary_file,
	"t:s"	=> \$taxon,
);

################

#if (not $blast2goannotfile or not $gff3file) {
if ((not $list_file) or (not $summary_file) or (not $taxon)) {
	print STDERR "\nThis script intersects output of scripts calculate_average_bootstrap_support.pl or calculate_average_branch_length.pl with SUMMARY file of RAXML tree construction to connect cluster ids with genes.

		Usage: link_tree_metrics_to_genes.pl -b <text_file_containing_cluster_ids_vs_metric> -s <tab_delimited_summary_file_of_trees> -t <abbreviation_of_taxon_in_question>\n\n";
	exit 1;
}


#print "gene_ID\tgene_length\tprotein_length\tavg_exon_length\tnumber_of_exons\tavg_intron_length\texon_proportion\tintron_proportion\n";
my $list_fh= &read_fh ($list_file);
my $summary_fh= &read_fh ($summary_file);
for (<$summary_fh>){
	chomp;
	push (@summary_array,$_);
}

while (<$list_fh>) {
	chomp;
	@list_array = split("\t");
#	print "$list_array[3]\n";
	for (@summary_array){
		@summary_array_split = split ("\t");
		if ($_ =~ $list_array[3]){
#			print $summary_array_split[4] . "\n";
			@taxa=split(",",$summary_array_split[4]);
			for (@taxa){
				if ($_ =~ $taxon){
#					print $list_array[0]."\t".$list_array[3]."\t".$_."\n";
					push(@out_array,$list_array[0]."\t".$list_array[3]."\t".$_);
				}
			}
			goto START;
		}
	}
	print "found NOTHING\n";
	START:
}

open (OUT,">$path/tree_metrics_vs_gene.tab") or die $!;
for (@out_array){
	print OUT $_."\n";
}
close OUT;
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

sub calc_GC {
	my $currentID=shift;
	my %GCcount=@_;
	my ($nucleotides_total,$GC_count)=(0,0);
#	print $currentID . "\n";
	foreach (keys %GCcount){
#		print "$_: $count{$_}\n";
		$nucleotides_total+=$GCcount{$_};
	}
#	print "total: $nucleotides_total\n";
	$GC_count += $GCcount{"G"} if ($GCcount{"G"});
	$GC_count += $GCcount{"C"} if ($GCcount{"C"});
#	print $GC_count."\n";
#	print "GC: ". $GC_count/$nucleotides_total . "\n";;
	return sprintf("%.3f",$GC_count/$nucleotides_total);
#	return $GC_count;
}
