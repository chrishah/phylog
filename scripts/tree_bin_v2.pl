#!/usr/bin/perl
#this script evaluates all treefiles (RAxML output Newick format) in a given directory (optionally a subset specified in a list file) and bins them into "good" - (i.e. no paralogs found), "good_paralog_mono" - (i.e. paralogous sequences are monophyletic) and "bad" datasets (i.e. paralogous sequences are paraphyletic).

use strict;
use warnings;
use POSIX qw(strftime);
use File::Copy qw(copy move);
use Getopt::Long;

my $fasta_dir = $ARGV[0];
my $input_filelist = $ARGV[1];
my @FILES;
my $exit;
my $minimal_alignment_length = 25;
#my $boostrap = 20;
my ($align, $trim, $phylo) = 0;
my @total; 
my @good_trees;
my @good_trees_paral_mono;
my @bad_trees_paral_para;
my $next_position;
my $prev_position;
my $size; 
my $neighbors;

my $USAGE = "\n	This script evaluates all treefiles (Newick format; filename *aln*) in a given directory (optionally a subset specified in a list file) and bins them into 
	\"good\" datasets (i.e. no paralogs found) 
	\"good_paralog_mono\" datasets (i.e. paralogous sequences are monophyletic) 
	\"bad\" datasets (i.e. paralogous sequences are paraphyletic)
		
		USAGE: ./tree_bin.pl <PATH/to/DIR> <optional-list-file-containing *.fasta files to be processed> \n\n";

if (!$ARGV[0]){
	print "$USAGE\n";
	exit;
}

opendir (DIR, "$fasta_dir") or die $!;
@FILES = readdir DIR;
closedir DIR;

if ($ARGV[1]){
	undef @FILES;
	open (INPUT, "<$input_filelist") or die $!;
	for (<INPUT>){
                push(@FILES, "$_");
        }
	close INPUT;

}
chdir "$fasta_dir" or die $!;

print strftime("%b %e %H:%M:%S", localtime) . "\n\n";
for (@FILES){
        chomp;
	my @IDs;
	my @ID;
	my %paralogs;
	my $treefile = $_;
	if ($_ ~~ /aln/){
       		print "\n-----\ncurrent tree is: $_\n\n";
	 	open (TREE,"<$_") or die $!;
		for (<TREE>) {
			chomp;
			my @treeID = split (/_/, $treefile);
			@treeID = splice @treeID, @treeID-2, $#treeID;
#			print "current @treeID\n";
			my @treeID_2 = split (/\./, $treeID[0]);
#			print "current: @treeID_2\n";
			unless (scalar @treeID_2 == 2){
				@treeID_2 = splice @treeID_2, -2, $#treeID_2;
			}
			my $treeID = join ".", @treeID_2;
			@treeID_2 = split (/\./, $treeID[1]);
#			print "current @treeID\n";
			$treeID = $treeID . "_" . $treeID_2[0];
			push (@total, $treeID);
#			print "\n-----\ncurrent tree: $treeID\n\n";
			print "$_\n\n";
			s/\(/ \( /g;
			s/\)/ \) /g;
			s/\,/ \, /g;
			s/:/ : /g;
			s/  / /g;
#			print "$_\n";
			my @elements = split / /;
			shift (@elements);
#			print "\n$elements[0]$elements[1]$elements[2]\n";
			for (@elements){
#				print "$_\n";
				if ($_ =~ /_/){
					@ID = split /_/;
					push (@IDs, $ID[0]);
				}
			}
#			print "IDs: " . scalar @IDs . "\n";
			for (@IDs){
#				print "$_\n";
			}
			my %count;
        		map { $count{$_}++} @IDs;
			my $total_size = scalar @IDs;
			foreach (keys %count){
 #               		print "$_ : $count{$_}\n";
				if ($count{$_} >= 2){
					$paralogs{$_} = $count{$_};
				}
			}
			if (!%paralogs){
				print "no paralogs found\n";
				push (@good_trees, $treeID);
			}else{
				print "paralogous sequences present in tree $treeID\n";
				foreach (keys %paralogs){
					my $current = $_;
					print "taxon $current is represented by $paralogs{$current} paralogs\n";
					$size = $paralogs{$current};
					$neighbors = 0;
					for (my $i = 0; $i < @IDs; $i++) {
						next unless $IDs[$i] ~~ /$current/;
#						print "$i: $IDs[$i]\n";
#						print "neighbors $neighbors\n";
#						print "current position $i\n";
						$next_position = $i + 1;
						$prev_position = $i - 1;
						if ($next_position == $total_size){
							$next_position = 0;
						}
						
#						print "previous $prev_position\n";
#						print "next $next_position\n";
						if ($IDs[$next_position] ~~ /$current/){
#							print "next position is a hit\n";
							$neighbors++;
						}else{
#							print "next position is no hit\n";	
						}
#						print "neighbors $neighbors\n";
						if ($IDs[$prev_position] ~~ /$current/){
#							print "previous position is a hit\n";
							$neighbors++;
						}else{
#							print "previous position is no hit\n";
						}
#						print "neighbors $neighbors\n";
					}
#					print "neighbors is: $neighbors\n";
#					my $sum = (2 + ($size - 2) * 2);
					my $sum = 2 * ($size - 1);
#					print "size is $size\n";
#					print "sum is $sum\n";
					if ($neighbors == (2 * ($size - 1))){
						print "taxon $current has passed the first monophyly test\n";
						my $parentheses = 0;
						my $count = 0;
						my $done = 0;
						my $number = $paralogs{$current};
#						print "number of paralogs: $number\n";
						for (my $i = 0; $i < @elements; $i++) {
#							print "$current count is $count\n";
#							print "current parentheses count is $parentheses\n";

							if (($count >= 1) && ($count < $number)){
#								print "arrived 1 < count < number: $elements[$i]\n";
								if ($elements[$i] =~ /\(/){
									$parentheses++;
								}elsif ($elements[$i] =~ /\)/){
									$parentheses--;
								}elsif ($elements[$i] =~ /$current/){
									$count++;		
#									print "found $count $current\n";
								}
							}elsif (($count == $number) && ($parentheses >= 1)){
								my $para = $parentheses;
								for (my $j = 1; $j <= $para; $j++){
#									print "cycle j is $j\n";
									$next_position = $i + ( 3 * $j - 1);
#									print "next position is: $next_position\n";
									if ($elements[$next_position] =~ /\)/){
#										print "parenthesis found in cycle j $j\n";
#										print "found parenthesis $j of $para\n";
										$parentheses--;
									}else {
										print "paraphyly alert 2\n";
										push (@bad_trees_paral_para, $treeID);
										goto END;
									}
								}
								if ($parentheses == 0){
									print "taxon $current has passed the second monophyly test\n\n";
									$done = 1;	
								}else {
									print "paraphyly alert 3\n";
									push (@bad_trees_paral_para, $treeID);
									goto END;
								}
							}elsif (($count == $number) && ($parentheses == 0)){
								print "taxon $current has passed the second monophyly test\n\n";
								$done = 1;	
							}

							last unless ($done == 0);
							next unless (($elements[$i] ~~ /$current/) && ($count == 0));
							$count++;
#							print "found $count $current\n";
							$prev_position = $i - 1;		
							if ($elements[$prev_position] =~ /\(/){
								$parentheses++;
							}else {
								print "paraphyly alert 1\n";
								push (@bad_trees_paral_para, $treeID);
								goto END;
							}
						}		
					}else{
						print "paraphyly detected for taxon $current in dataset $treeID\n";
						push (@bad_trees_paral_para, $treeID);
						goto END;
					}
					
				}
				push (@good_trees_paral_mono, $treeID);
				print "all paralogs in dataset $treeID seem to be monophyletic\n";
				END:
			}
		}
		close TREE;
	}
}
unless (!@good_trees){
	open (GOOD,">good_IDs.list") or die $!;
	for (@good_trees){
		print GOOD "$_\n";
	}
	close GOOD;
}
unless (!@good_trees_paral_mono){
	open (GOOD_PARALOG,">good_but_paralogs.list");
	for (@good_trees_paral_mono){
		print GOOD_PARALOG "$_\n";
	}
close GOOD_PARALOG;
}
unless (!@bad_trees_paral_para){
	open (BAD,">bad_IDs.list") or die $!;
	for (@bad_trees_paral_para){
		print BAD "$_\n";
	}
	close BAD;
}
my $percent_good = (scalar @good_trees) / scalar @total;
my $percent_good_paral_mono = (scalar @good_trees_paral_mono) / scalar @total;
my $percent_bad_paral_para = (scalar @bad_trees_paral_para) / scalar @total;  
print 	"\n\n	-------- SUMMARY --------
	total number of clusters: " . scalar @total . "
	clusters without paralogs: " . scalar @good_trees . " (" . sprintf ("%.2f", (scalar @good_trees / scalar @total) * 100) . " %)
	clusters with monophyletic paralogs: " . scalar @good_trees_paral_mono . " (" . sprintf ("%.4f", (scalar @good_trees_paral_mono / scalar @total)) * 100 . " %)
	clusters with paraphyletic paralogs: " . scalar @bad_trees_paral_para . " (" . sprintf ("%.2f", (scalar @bad_trees_paral_para / scalar @total) * 100) . " %)\n\n";
