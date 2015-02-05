#!/usr/bin/perl

use strict;
use warnings;
use POSIX qw(strftime);
use File::Copy qw(copy move);
use Getopt::Long;


############################################
#set path to external scripts/programs first
############################################
my $clustalo="clustalo-1.1.0-linux-64";
my $ALICUT="ALICUT_V2.3.pl";
my $ALISCORE="Aliscore.02.2.pl";
my $select_contigs="select_contigs.pl";
my $fasta2phylip="fasta2phylip.pl";
my $RAxML="raxmlHPC-PTHREADS -T 2";
my $Protein_model_select="ProteinModelSelection.pl";
if ((!$clustalo)||(!$ALICUT)||(!$ALISCORE)||(!$select_contigs)||(!$fasta2phylip)||(!$RAxML)||(!$Protein_model_select)){
	print "\n\nplease set paths to external scripts/programs first\n\n";
	exit;
}
#################################################

my $fasta_dir = $ARGV[0];
my $input_filelist = $ARGV[1];
my (@stats, @FILES, @taxa, @IDs, @full_fasta, @full_taxa, @paralogIDs, @paralog_taxa, @orthologIDs, @alternative_alignments, @alignment_length);
my ($summary, $taxa_sum, $exit, $alignment_length, $unique_taxa, $model);
my %count;
my $minimal_alignment_length = 50;
my $bootstrap = 20;
my ($help, $align, $mono_mask, $trim, $notree, $phylo, $current) = (0, 0, 0, 0, 0, 0, 0);
my $USAGE = "\nThis script is a wrapper for phylogenetic inference. If no options are specified it will perform (for every *.fasta file in the specified directory):

	-protein sequence alignment using clustalo
	-alignment trimming using ALISCORE and ALICUT
	-monophyly masking, i.e. if more than one sequence per taxon is present it will choose one representative as described in Hahn et al. 2014
	-find best fitting model of protein evolution for the remaining alignment
	-infer ML tree using the infered model using RAxML

Each of these steps can be individually omitted via the options.


	USAGE: ./process_genes.pl <PATH/to/DIR> <options>

	options:\n
	--help		show this infomation
	--noalign	omit alignment
	--nomask	omit monophyly masking
	--notrim	omit alignment trimming
	--notree	omit tree inference (but do model prediction)
	--nophylo	omit phylogenetic analysis (i.e. model prediction and tree inference)
	--bootstrap	number of bootstrap replicates to perform in phylogenetic analysis (default=20)\n";
		

GetOptions (	"help!" => \$help,
		"noalign!" => \$align,
		"nomask!" => \$mono_mask,
                "notrim!" => \$trim,
		"notree!" => \$notree,
		"nophylo!" => \$phylo,
		"bootstrap:i" => \$bootstrap) or die "Incorrect usage!\n$USAGE";


if ((!$ARGV[0]) || ($help)){
	print "$USAGE\n";
	exit;
}

print "\nparameters (0=off; 1=on):
	align: $align
	trim: $trim
	monomask: $mono_mask
	notree: $notree
	nophylo: $phylo\n";

opendir (DIR, "$fasta_dir") or die $!;
@FILES = readdir DIR;
closedir DIR;

chdir "$fasta_dir" or die $!;
$summary = "\#cluster name    number of taxa  alignment length        model   list of taxa";
push (@stats, $summary); 
print strftime("%b %e %H:%M:%S", localtime) . "\n\n";
for (@FILES){
        chomp;
        if ($_ =~ /fasta$/){
                print "------\ncurrent: $_\n";
                my @file = split (/\./, $_);
		$summary = $_;
#               print "$file[0]\n";
#               print "$file[1]\n";
#               print $prefix . "\n";
		my $prefix = "$file[0]" . "." . "$file[1]";
		mkdir "$prefix\_processed" or die $!;
		open (FASTAFILE,"<$_") or die $!;
		open (FASTAFILEWITHOUTSTAR,">$prefix\_processed/$_") or die $!;
		for(<FASTAFILE>){
			chomp;
			if ($_ =~ /^>/){
				s/\|/_/g;
				print FASTAFILEWITHOUTSTAR $_ . "\n";
			}
			else{
				s/\*//g;
				s/U/C/g;
				print FASTAFILEWITHOUTSTAR $_ . "\n";
			}
		}
		close FASTAFILEWITHOUTSTAR;
		#print $prefix . "\n";
		chdir "$prefix\_processed" or die $!;	
		if (!$align){
			print "running clustal omega\n";
			my @clustal_output = `$clustalo -i $prefix.fasta -o $prefix.aln.fasta`;
			$exit = $? >> 8;
        		unless ($exit == 0){
                		print "\nsome kind of ERROR occured while running clustalo\n";
        			goto END;
			}
			print strftime("%b %e %H:%M:%S", localtime) . "\n\n";
                }else {
			print "alignment omitted\n";
		}
		$prefix = "$file[0]" . "." . "$file[1]" . ".aln";
		
		if (!$mono_mask){
			open (BASEFASTA,"<$prefix.fasta") or die $!;
			undef @full_fasta;
 			undef @full_taxa;
			undef @taxa;
			undef @IDs;
			my %taxa_hash = ();
			for (<BASEFASTA>){
				chomp;
				push (@full_fasta, $_);
				if ($_ =~ /^>/){
					s/^>//g;
					@taxa = split (/_/, $_);
					$taxa_hash{$taxa[0]} = $taxa[1];
					push (@IDs, $taxa[0]);
					push (@full_taxa, $_);
				}
			}
#			foreach (keys %taxa_hash){
#				print "$_\n";
			close BASEFASTA;
			%count = ();
			map { $count{$_}++} @IDs;
			my $size = scalar keys %count;
			print "number of taxa: $size\n";
			foreach (keys %count){
				$current = $_;
				if ($count{$_} >= 2){
					print "$_ : $count{$_}\n";
					push (@paralog_taxa, $_);
					for (@full_taxa){
#						print "compare $current to $_\n";
						if ($_ =~ /^$current/){
							push (@paralogIDs, $_);
						}
					}	
				}else {
					print "$_ : $count{$_}\n";
					for (@full_taxa){
						if ($_ =~ /^$current/){
							push (@orthologIDs, $_);
						}
					}
				}
			}
			my $number_of_orthologs = scalar @orthologIDs;
			print "number of unique taxa: $number_of_orthologs\n";
#			print "orthologs:\n";
#			for (@orthologIDs){
#				print "$_\n";
#			}
#			print "paralogs\n";
			undef @alternative_alignments;
			for (my $i = 0; $i < @paralogIDs; $i++){
#				print "push $paralogIDs[$i] to array\n";
				push (@orthologIDs, $paralogIDs[$i]);
				for (@orthologIDs) {
					open (LIST,">>orthologs_plus_$paralogIDs[$i]");
#					print "$_\n";
					print LIST "$_\n";
				}		
				`$select_contigs -n orthologs_plus_$paralogIDs[$i] $prefix.fasta orthologs_plus_$paralogIDs[$i].aln.fasta`;
				unlink "orthologs_plus_$paralogIDs[$i]";
				close LIST;
				push (@alternative_alignments, "orthologs_plus_$paralogIDs[$i].aln.fasta");
				pop (@orthologIDs);	
			}
			print "alternative alignments: " . scalar @alternative_alignments . "\n";
#			print "trim ambiguously aligned positions from alignment\n";
			for (@alternative_alignments){
				&ALISCORE($_);	
			}
			&ALICUT;
#			}
#		}else {
#			print "monophyly masking omitted\n";
#		}

			print "choosing paralogs\n";
#			foreach (keys %count){
			for (@paralog_taxa){
				$current = $_;
				my %hash = ();
				my $longest_alignment_length = 0;
				my $longest_alignment_ID = 0;
#				print "current Taxon: $current\n";
				for (@alternative_alignments){
					my @tmp	= split (/\./, $_);
					pop (@tmp);
					$prefix = join (".", @tmp);	
					if ($_ ~~ /$current/){
						($unique_taxa, $alignment_length, $taxa_sum) = &convert_phy_length("ALICUT_$prefix");
	
						if ($alignment_length >= $longest_alignment_length){
							@tmp = split (/\./, $_);
							@tmp = splice @tmp, 0, -2;
							my $tmp = join ('.', @tmp);
							@tmp = split (/_/, $tmp);
							@tmp = splice @tmp, 2, $#tmp;
							$tmp = join ('_', @tmp);
							$longest_alignment_ID = $tmp;
							$longest_alignment_length = $alignment_length;
						}
					}
				}
				print "The longest alignment is $longest_alignment_ID with length $longest_alignment_length\n";
				push (@orthologIDs, $longest_alignment_ID);	
				open (LONG,">best_alignment.list");
				for (@orthologIDs){
#					print "$_\n";
					print LONG "$_\n";
				}
				close LONG;
				$prefix = "$file[0]" . "." . "$file[1]" . ".aln";	
			}
			
			`rm ./*orthologs_*`;
			`$select_contigs -n best_alignment.list $prefix.fasta $prefix.opt.fasta`;
			$prefix = $prefix . ".opt";
		}else {
			print "monophyly masking omitted\n";
		}

		if (!$trim){
#			&ALISCORE("$prefix.opt.fasta");
			&ALISCORE("$prefix.fasta");
			&ALICUT;
#			$prefix = "ALICUT_" . "$file[0]" . "." . "$file[1]" . ".aln.opt";		
			$prefix = "ALICUT_" . "$prefix";		
			($unique_taxa, $alignment_length, $taxa_sum) = &convert_phy_length("$prefix");
			if ($alignment_length < $minimal_alignment_length){
				print "the alignment is shorter than $minimal_alignment_length and will not be further processed\n";
				move ("$prefix.phy", "$prefix.phy.tooshort") or die $!;
				$model = "-";
				goto END;
			}
		}else {
			print "alignment trimming omitted\n";
		}
		
		if (!$phylo){
			if (($trim)&&($mono_mask)){
				print "notrim && nomonomask\n";
				
				$prefix = substr($summary,0,-6);
				print "prefix: $prefix\n";
				($unique_taxa, $alignment_length, $taxa_sum) = &convert_phy_length("$prefix");
				if ($alignment_length < $minimal_alignment_length){
	                                print "the alignment is shorter than $minimal_alignment_length and will not be further processed\n";
        	                        move ("$prefix.phy", "$prefix.phy.tooshort") or die $!;
                	                $model = "-";
                        	        goto END;
				}
                        }
			$model = &RAxML_model($prefix);

#			}
			if (!$notree){
				&RAxML_tree($model, $bootstrap, $prefix);
			}else{
				print "best fitting AA model was found to be $model\n";
				print "ML tree inference omitted\n";
			}
			
		}else {
			print "phylogenetic analyses omitted\n";
			$model = "-";
		}
		
		END:
		undef @orthologIDs;
		undef @paralog_taxa;
		undef @paralogIDs; 
		chdir "../" or die $!;
		print strftime("%b %e %H:%M:%S", localtime) . "\n\n";

		$summary = $summary . "\t" . $unique_taxa . "\t" . $alignment_length . "\t" . $model . "\t" . $taxa_sum;
		push (@stats, $summary);
	}
}
open (STATS,">stats.tab") or die $!;
map { print STATS "$_\n" } @stats;
close STATS;


sub ALICUT{
        print "\nrun Alicut\n\n";
        my @Alicut_output = `$ALICUT -s`;
        $exit = $? >> 8;
        unless ($exit == 0){
                print "\nsome kind of ERROR occured while running Alicut\n";
                goto END;
        }
        open (OUTPUT,">Alicut.log") or die $!;
        for (@Alicut_output){
                print OUTPUT $_ . "\n";
        }
        close OUTPUT;
        print strftime("%b %e %H:%M:%S", localtime) . "\n\n";
}

sub ALISCORE{
        print "run Aliscore on $_[0]\n";
        my @Aliscore_output = `$ALISCORE -N -r 200000000000000000 -i $_[0]`;
        $exit = $? >> 8;
        unless ($exit == 0){
                print "\nsome kind of ERROR occured while running Aliscore\n";
        }
        open (OUTPUT,">Aliscore-$_[0].log") or die $!;
        for (@Aliscore_output){
                print OUTPUT $_[0] . "\n";
        }
        close OUTPUT;
}

sub convert_phy_length{
        print "convert $_[0] to *.phylip\n";
        my @conversion = `$fasta2phylip $_[0].fasta $_[0].phy`; #convert *.fasta to *.phylip
        $exit = $? >> 8;
        unless ($exit == 0){
                print "\nsome kind of ERROR occured while running fasta2phylip conversion\n";
                goto END;
        }
	my (@phylip, @ID, @taxon, @new_phylip, @alignment_length) = ();
        open (PHYLIP,"<$_[0].phy") or die $!;
	map {chomp; push (@phylip, $_)} <PHYLIP>;
	close PHYLIP;
	my $bad_count = 0;
	for (my $i = 0; $i < @phylip; $i++){
		if ($i == 0){
			@alignment_length = split (/ /, $phylip[$i]);
			chomp @alignment_length;
 #       		print "number $alignment_length[0] length $alignment_length[1]\n";
		}elsif ($i >= 1){
			@ID = split (/\t/, $phylip[$i]); 
#			print "current line: $phylip[$i]\n";
			if ($ID[1] =~ /[A-Z]/){
#				print "ok\n";
				push (@new_phylip, $phylip[$i]);
				push (@taxon, $ID[0]); 
			}else{
				$bad_count++;	
#				print "not ok. count: $bad_count\n";
			}
		}
	}
	if ($bad_count >= 1){
		open (PHYLIP,">$_[0].phy") or die $!;
#		my $length = $alignment_length[0];
#		print "original count = $alignment_length[0]\n";
		$alignment_length[0] -= $bad_count;
#		print "new count = $alignment_length[0]\n";
		print PHYLIP "$alignment_length[0] $alignment_length[1]\n";
		map {chomp; print PHYLIP "$_\n";}@new_phylip;
		close PHYLIP;
	}
#	my @alignment_length = split (/ /, <PHYLIP>);
#	chomp @alignment_length;
#	my (@sequence, @taxon) = ();
#	map {my @ID = split /\t/; push (@taxon, $ID[0]); push (@sequence, $ID[1])} <PHYLIP>;
        print "alignment length after trimming is: $alignment_length[1]\n";
	@taxon = sort @taxon;
	my $taxa = join (",", @taxon);
        return ($alignment_length[0], $alignment_length[1], $taxa);
}

sub RAxML_model{
	print "find best fitting protein model using ProteinModelSelection.pl script\n";
	my @RAxML_AA_Test_output = `$Protein_model_select $_[0].phy`; #find best fitting model
	$exit = $? >> 8;
	unless ($exit == 0){
		print "\nsome kind of ERROR occured while running RAxML_AA_Test.pl script\n";
		goto END;
	}

	open (MODEL,">$_[0].bestmodel");
	print MODEL @RAxML_AA_Test_output;
	close MODEL;
#	print "@RAxML_AA_Test_output\n";
	print strftime("%b %e %H:%M:%S", localtime) . "\n\n";
	my @model;
	for (@RAxML_AA_Test_output){
		chomp;
		if ($_ =~ /Model/){
			@model = split (/: /, $_);
		}
	}
	return $model[1];
}

sub RAxML_tree{
	print "running RAxML-HPC using the determined model: ($_[0]) and $_[1] bootstrap replicates\n";
        my $RAxMLmodel = "PROTGAMMA" . $_[0];
        my @RAxML_run = `$RAxML -f a -m $RAxMLmodel -p 12345 -x 12345 -# $_[1] -s $_[2].phy -n $_[2]`;
        unless ($exit == 0){
                print "\nsome kind of ERROR occured while running RAxML\n";
                goto END;
        }

}
