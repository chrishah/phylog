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
	print STDERR "\nThis script intersects a list of gene models from a given taxon with the orthomcl output file (groups.txt) and outputs all gene IDs for a specified taxon that are present in the same clusters.

		Usage: find_taxon_gene_models_overall_sharing_cluster.pl -b <list of gene models of a taxon> -s <groups.txt> -t <abbreviation_of_taxon_in_question>\n\n";
	exit 1;
}


#print "gene_ID\tgene_length\tprotein_length\tavg_exon_length\tnumber_of_exons\tavg_intron_length\texon_proportion\tintron_proportion\n";
my $list_fh= &read_fh ($list_file);
my $summary_fh= &read_fh ($summary_file);
my (@taxon_ids,@search_taxon_all_ids,@cluster);
my $i;
for (<$summary_fh>){
	chomp;
	push (@summary_array,$_);
}

for (<$list_fh>){
	chomp;
	push(@list_array,$_);
}
my $search_taxon = substr $list_array[1], 0,3;
#print "$search_taxon\n";
for (@list_array) {
	chomp;
	my @search_taxon_id=split("_");
	print "search for: $search_taxon_id[0]\|$search_taxon_id[1]\n";
#	for (@summary_array){
#		chomp;
#		my $line = $_." ";
	for ($i=0; $i<@summary_array; $i++){
                chomp;
                my $line = $summary_array[$i]." ";
#		print "$line\n";
#		if ($line =~ /$search_taxon_id /){
		if ($line =~ /$search_taxon_id[0]\|$search_taxon_id[1] /){
			@summary_array_split = split (": ",$summary_array[$i]);
			push(@cluster,$summary_array_split[0]);
			print "found in cluster: $summary_array_split[0]\n";
#			print "$line\n";
			my @IDs=split(" ",$summary_array_split[1]);
#			@taxa=split(",",$summary_array_split[4]);
#			print "@taxa\n";
			for (@IDs){
				if ($_ =~ $taxon){
					print "found: $_\n";
#					my @cluster_temp = split ('\.',$summary_array_split[0]);
#					my @cluster = splice @cluster_temp, 0, 2;
#					my $cluster = join(".",@cluster);
#					my $cluster = $summary_array_split[0];
#					print "$cluster\n";
#					push(@out_array,$summary_array_split[0]."\t".$search_taxon_id[0]."\|".$search_taxon_id[1]."\t".$_);
					push(@taxon_ids,$summary_array_split[0]."\t".$_);
				}elsif($_ =~ $search_taxon_id[0]){
					print "found: $_\n";
					push(@search_taxon_all_ids,$summary_array_split[0]."\t".$_);
				}
			}
			last;
		}

	}
	splice(@summary_array,$i,1);
	print "array contains ".scalar @summary_array." elements\n";
}

open (OUT_1,">$path/$search_taxon.list") or die $!;
open (OUT_2,">$path/$taxon.list") or die $!;
open (OUT_clusters,">$path/clusters.list") or die $!;

#print "before ". scalar @search_taxon_all_ids."\n";

#my @uniq_search_taxon = &remove_redundant(@search_taxon_all_ids);
#my @uniq_taxon = &remove_redundant(@taxon_ids);
#my @uniq_clusters = &remove_redundant(@cluster);

#print "after ". scalar @uniq_search_taxon."\n";

for (@search_taxon_all_ids){
	print OUT_1 $_."\n";
}
for (@taxon_ids){
	print OUT_2 $_."\n";
}
for (@cluster){
	print OUT_clusters $_."\n";
}


################################

sub remove_redundant{
	my %seen;
	$seen{$_}++ for @_;
	my @unique = keys %seen;
	return @unique;
}
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
