#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);

my ($clusterfile, $listfile) = ("","");
my (@groups, @array);
#my ($blast2goannotfile, $gff3file, $direction, $gene_length, $id, $exon_length, $intron) =("","","","","","","");
#my ($check, $exon_length_cummulative, $intron_length_cummulative, $average_intron) = 0;
#my (@exon_per_gene, @intron_per_gene, @exons_total, @average_exon_length, @average_exon_per_gene, @average_intron_length, @total_exon_length, @total_intron_length);
GetOptions (
#	"blast2goannotfile:s" => \$blast2goannotfile,
	"clusterfile:s"    => \$clusterfile,
	"listfile:s"	=> \$listfile,
);

################

#if (not $blast2goannotfile or not $gff3file) {
if (not $clusterfile or not $listfile){
	print STDERR "\nThis script gets sequence IDs from clusters. It needs to be fed with a list of clusters and a file containing a list of sequences per clusters in the format as produced by mcl, i.e. groups.txt

		Usage: get_sequence_IDs_from_cluster_IDs.pl -c <groups.txt> -l <list>\n\n";
	exit 1;
}

#my $blast2goannotfh = &read_fh ($blast2goannotfile);

#my %note;
#my %go;
#my %ec;
#
#while (<$blast2goannotfh>) {
#	chomp;
#	my @F = split /\t/;
#	$note{$F[0]}.=$F[2]   if exists $F[2];
#	$go  {$F[0]}{$F[1]}=1 if $F[1]=~/^GO:/;
#	$ec  {$F[0]}{$F[1]}=1 if $F[1]=~/^EC:/;
#}

################

my $listfh = &read_fh ($listfile);
my $clusterfh = &read_fh ($clusterfile);

for (<$clusterfh>){
	chomp;
	push (@groups, $_);
}

while (<$listfh>) {
	chomp;
	my $current = $_;
	for (@groups){
		chomp;
		if ($_ =~ /$current/){
			undef @array;
			@array = split / /;
			shift @array;
			for (@array){
				chomp;
				print "$_\n";
			}
		}
	}
}

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
