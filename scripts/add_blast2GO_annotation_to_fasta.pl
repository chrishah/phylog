#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);
use Cwd;

my $path=getcwd;
my $fasta_file;
my $annotation_file;
my $out_file;
my (@fasta_array, @annotation_array, @out_array);
my (@temp_1, @temp_2, @temp_3);

GetOptions (
#	"blast2goannotfile:s" => \$blast2goannotfile,
	"f:s"    => \$fasta_file,
	"a:s"	=> \$annotation_file,
	"o:s"	=> \$out_file,
);

################

#if (not $blast2goannotfile or not $gff3file) {
if ((not $fasta_file) or (not $annotation_file) or (not $out_file)) {
	print STDERR "\nThis script replaces the header in a multifasta file (header format: somtething|ID) with annotations produced by blast2go (contained in a fasta like list, with each line looking like >something|ID|annoation). 
		Usage: add_blast2GO_annotation_to_fasta.pl -f multifasta.fasta -a annotations_list -o out.fasta\n\n";
	exit 1;
}


#print "gene_ID\tgene_length\tprotein_length\tavg_exon_length\tnumber_of_exons\tavg_intron_length\texon_proportion\tintron_proportion\n";
my $fasta_fh= &read_fh ($fasta_file);
my $annotation_fh = &read_fh ($annotation_file);
for (<$fasta_fh>){
	chomp;
	push (@fasta_array,$_);
}

for (<$annotation_fh>){
	chomp;
	if ($_ =~ /^>/){
		push(@annotation_array, $_);
	}
}

for (@fasta_array) {
	chomp;
	if ($_ =~ /^>/){
		my $test = 0;
#		print "fasta file: $_\n";
		@temp_1 = split (/\|/);
		for (my $i=0; $i<@annotation_array; $i++){
			@temp_2 = split (/\|/, $annotation_array[$i]);
			if ($temp_2[1] eq $temp_1[1]){
#				print "current ID: $temp_1[0]|$temp_1[1]\n";
#				print "array size before: ".scalar @annotation_array."\n";
				$test = 1;
#				print "annotation: $temp_2[1]\n";
#				print "$_\n";
				push (@out_array,$annotation_array[$i]);
				my $spliced = splice @annotation_array, $i, 1;
#				print "spliced: $spliced\n";
#				print "array size after: ".scalar @annotation_array."\n";
				last;
			}
		}	
		if (!$test){
			print "somethings wrong with $_\n";
			exit;
		}
	}else{
#		print "$_\n";
		push (@out_array,$_);
	}
	END:
}

open (OUT,">$out_file") or die $!;
for (@out_array){
	print OUT "$_\n";
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

