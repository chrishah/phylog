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
	print STDERR "\nThis script intersects a list of gene models from a given taxon with Tree-stats summary file and outputs all gene IDs for a specified taxon that are present in these clusters.

		Usage: find_taxon_gene_models_sharing_cluster.pl -b <list of gene models of a taxon> -s <tab_delimited_summary_file_of_trees> -t <abbreviation_of_taxon_in_question>\n\n";
	exit 1;
}


#print "gene_ID\tgene_length\tprotein_length\tavg_exon_length\tnumber_of_exons\tavg_intron_length\texon_proportion\tintron_proportion\n";
my $count=0;
my $list_fh= &read_fh ($list_file);
my $summary_fh= &read_fh ($summary_file);
for (<$summary_fh>){
	chomp;
	push (@summary_array,$_);
}

for (<$list_fh>){
	chomp;
	push(@list_array,$_);
}
for (@list_array) {
	chomp;
	my $search_taxon_id=$_;
	print "$search_taxon_id\n";
	for (@summary_array){
		if ($_ =~ /$search_taxon_id,/){
			@summary_array_split = split ("\t");
			@taxa=split(",",$summary_array_split[4]);
#			print "@taxa\n";
			for (@taxa){
				if ($_ =~ $taxon){
#					print $_."\n";
					my @cluster_temp = split ('\.',$summary_array_split[0]);
					my @cluster = splice @cluster_temp, 0, 2;
					my $cluster = join(".",@cluster);
#					print "$cluster\n";
					push(@out_array,$cluster."\t".$search_taxon_id."\t".$_);
					$count++;
				}
			}
		}
#		if ($search_taxon_id =~ /$_/){
#			for (@taxa){
#				if ($_ =~ $taxon){
#					my @cluster_temp = split (".",$summary_array_split[0]);
#					my @cluster = splice @cluster_temp, 0, 2;
#					my $cluster = join(".",@cluster);
#					push(@out_array,$cluster."\t".$search_taxon_id."\t".$_);
#				}
#			}
#			goto START;
#		}
	}
#	print "found NOTHING\n";
#	START:
}

open (OUT,">$path/OrthoClusterID_shared_genes.tab") or die $!;
for (@out_array){
	print OUT $_."\n";
}
close OUT;
print "found corresponding $taxon id in $count clusters\n";

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
