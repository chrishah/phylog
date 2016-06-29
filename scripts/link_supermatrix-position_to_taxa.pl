#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);
use Cwd;

my $path=getcwd;
my $partition_file;
my $summary_file;
my $annotation_file;
my $out;
my $taxon;
my $annot;
my (@partition_array, @summary_array, @summary_array_split,@annotation_array,@out_array,@taxa);
my (@temp_1, @temp_2, @temp_3);

GetOptions (
#	"blast2goannotfile:s" => \$blast2goannotfile,
	"p:s"    => \$partition_file,
	"s:s"	=> \$summary_file,
	"o:s"	=> \$out,
	"a:s"	=> \$annotation_file,
	"t:s"	=> \$taxon,
);

################

#if (not $blast2goannotfile or not $gff3file) {
if ((not $partition_file) or (not $summary_file) or (not $out) or (not $annotation_file) or (not $taxon)) {
	print STDERR "\nThis script intersects a partition.file (format as required for RAxML, format: MODEL, cluster_name = start_pos-end_pos), a Tree-stats summary file and outputs a file (name specified by the user) of format: column1= position_from => to; column2= # of taxa; column3= comma separated list of taxaIDs, and a annotated fasta file in the format >bla|bla|annotation. Also required is to specify a taxon id (3-letter, e.g. Gsa), that the scripts needs to link the info together.

		Usage: link_supermatrix-position_to_taxa.pl -p <part.raxml> -s <tab_delimited_summary_file_of_trees> -a annotated.fasta -t Gsa -o output.file\n\n";
	exit 1;
}


#print "gene_ID\tgene_length\tprotein_length\tavg_exon_length\tnumber_of_exons\tavg_intron_length\texon_proportion\tintron_proportion\n";
my $partition_fh= &read_fh ($partition_file);
my $summary_fh= &read_fh ($summary_file);
my $annotation_fh = &read_fh ($annotation_file);
for (<$summary_fh>){
	chomp;
	push (@summary_array,$_);
}

for (<$partition_fh>){
	chomp;
	push(@partition_array,$_);
}
for (<$annotation_fh>){
	chomp;
	if ($_ =~ /^>/){
		push(@annotation_array, $_);
	}
}

for (@partition_array) {
	chomp;
	@temp_1 = split (" ");
	my $cluster_ID = $temp_1[1].".";
	$temp_1[3] =~ s/-/ => /g;	
	

#	my @temp = split ("\t");
#	print "$temp[0]\n";
#	my @temp2 = split ("_", $temp[0]);
#	print "$temp2[1]\t$temp2[2]\n";
#	my $cluster_name = join ("_", $temp2[1], $temp2[2]);
#	print "$cluster_name\n";
#	my @temp3 = split (/\./, $cluster_name);
#	print "$temp3[0]\t$temp3[1]\n";
#	my $cluster_ID = join (".", $temp3[0], $temp3[1]);
#	$cluster_ID .= ".";
#	print "$cluster_ID\n";
#	goto END;
	for (@summary_array){
		if ($_ =~ /$cluster_ID/){
			@summary_array_split = split ("\t");
			my @gene_model_id = split (",",$summary_array_split[4]);
#			print $temp_1[3]."\t".$summary_array_split[2]."\t".$summary_array_split[1]."\t".$summary_array_split[4]."\n";			
			for (@gene_model_id){
				if ($_ =~ /$taxon/){
#					print "1st: $_\n";
					@temp_2=split("_");
#					$temp_2[1] = "|".$temp_2[1]."|";	
#					print "id is: $temp_2[1]\n";
					for (@annotation_array){
						@temp_3 = split (/\|/);
#						print "id: $temp_3[1]\n";
#						if (($temp_3[0] eq $taxon) and ($temp_3[1] == $temp_2[1])){
						if (($temp_3[1] == $temp_2[1])){
#						print "annotation: $_\n";
#						if ($_ =~ /$temp_2[1]/){
#							print "2nd: $_\n";
#							@temp_4 = split ("|");
							$annot = $temp_3[2];
							if ((!$annot) or ($annot =~ /NA/)){
								print "could not find or having problems with gene model $_\n";
								exit;
							}
						}
					}
				}
			}
#			push (@out_array, $temp[1]."\t".$summary_array_split[2]."\t".$summary_array_split[1]."\t".$summary_array_split[4]);
			push (@out_array, $temp_1[3]."\t".$summary_array_split[2]."\t".$summary_array_split[1]."\t".$summary_array_split[4]."\t". $annot);
#			print $temp[1]."\t".$summary_array_split[2]."\t".$summary_array_split[1]."\t".$summary_array_split[4]."\n";

			goto END;
		}
	}
#	print "found NOTHING\n";
	push (@out_array, $cluster_ID."\tnothing found");
	END:
}

open (OUT,">$path/$out") or die $!;
print OUT "#position in matrix\tlength\tnumber of taxa\tlist of taxon IDs\tannotation of $taxon gene model\n";
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

