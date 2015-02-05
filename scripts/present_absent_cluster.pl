#!/usr/bin/env perl

use strict;
use warnings;

if ((!$ARGV[0]) or (!$ARGV[1]) or (!$ARGV[2]) or (!$ARGV[3])){
	print STDERR "\nThis script checks a groups.txt OrhtoMCL output and finds the clusters that are presnet/absent in all taxa specified by the user.

                Usage: present_absent_cluster.pl <groups.txt> <list containing taxa that should be contined in the cluster> <list of taxa that should now be contined> <output file name>

		example: present_absent_cluster.pl groups.txt.gz present_example.txt absent_example.txt output.txt
\n\n";
        exit 1;
}


my $groups_fh=&read_fh ($ARGV[0]);
my $contained_list_fh=&read_fh ($ARGV[1]);
my $absent_list_fh=&read_fh ($ARGV[2]);
my $output=$ARGV[3];
#my $taxon_id="Gsa";
my (@groups_array,@contained_list_array,@absent_list_array);

while(<$groups_fh>){
	chomp;
	push(@groups_array,$_);
}
while(<$contained_list_fh>){
	chomp;
	push(@contained_list_array,$_);
}
while(<$absent_list_fh>){
	chomp;
	push(@absent_list_array,$_);
}

my %count;
my @cluster_array;
my $groups_line;

open (OUT,">$output") or die $!;
#open (OUT_length,">tree_metric_vs_xon_length.csv") or die $!;

for(@groups_array){
	$groups_line=$_;
#	my @cluster=split(": ");
#	my $cluster=shift(@cluster);
#	print "cluster name is: $cluster\n";
	my @gene_id=split(" ");
	my $cluster=shift(@gene_id);
	$cluster = substr $cluster,0,-1;
#	print "\nstart cluster $cluster\n";
#	for (@gene_id){
#		print "$_\n";
#	}
	undef %count;
	map { s/\|.*//g; $count{$_}++} @gene_id;	
#	foreach (keys %count){
#		print "$_: $count{$_}\n";
#	}
	my $j = 0;
	for (my $i=0; $i<@contained_list_array; $i++){
#		print "actual taxon id is $contained_list_array[$i]\n";
#		print "number is $count{Hsa}\n";
		if (!$count{$contained_list_array[$i]}){
#			print "cluster $cluster does not contain taxon $contained_list_array[$i]\n";
			last;
		}else{
#			print "$contained_list_array[$i] contained in cluster!\n";
			$j++;
		}
		if ($j==scalar @contained_list_array){	
#			print "all mandatory taxa contained! checking for absence..\n";
			for (@absent_list_array){
				if ($count{$_}){
#					print "cluster $cluster contains $_\n";
					goto JUMP;
				}
			}
			push (@cluster_array, $groups_line);
#			print "cluster $cluster seems to be ok!\n";
			foreach (keys %count){
#				print "$_: $count{$_}\n";
			}
		}
		JUMP:
	}
}

for (@cluster_array){
	print OUT "$_\n";
}
print "total number of clusters meeting the criteria: ". scalar @cluster_array. "\n";

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
