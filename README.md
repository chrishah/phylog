Introduction
------------

This repository contains useful script for the preparation of data for phylogenomic analyses.

The scripts in this repository were used in the analyses presented in [Hahn et al. 2014](http://gbe.oxfordjournals.org/content/early/2014/04/13/gbe.evu078.short?rss=1 "Hahn et al. 2014 at GBE"). They were tailored specifically for this study, but might be useful for other users as well. Obviously, I won't guarantee that they will work for your data. If you find them useful, I would appreciate if you would cite our [paper](http://gbe.oxfordjournals.org/content/early/2014/04/13/gbe.evu078.short?rss=1 "Hahn et al. 2014 at GBE"). 

If you have questions/comments don't hesitate to contact me at: c.hahn@hull.ac.uk

A step by step guide through some parts of the process described in [Hahn et al. 2014](http://gbe.oxfordjournals.org/content/early/2014/04/13/gbe.evu078.short?rss=1 "Hahn et al. 2014 at GBE") is given below.


Prerequisites
-------------

A number of external scripts/programs are expected to be in your path:

- sort_contigs.pl - script available from [here](http://www.genome.ou.edu/informatics.html) (required by sort_clusters_v2.pl if you want to extract the sequences into a fasta file - see below)
- fasta2phylip.pl - script available from [here](https://github.com/chinchliff/physcripts/blob/master/fasta2phylip.pl) (required by process_genes.pl)
- clustal omega - available from [here](http://www.clustal.org/omega/) (required by process_genes.pl) - scripts were tested with clustal omega 1.1 
- ALICUT - available from [here](http://zfmk.de/web/Forschung/Abteilungen/AG_Wgele/Software/Utilities/index.en.html)(required by process_genes.pl) - scripts were tested with ALICUT_v2.3
- ALICSCORE - availalbe from [here](http://zfmk.de/web/Forschung/Abteilungen/AG_Wgele/Software/Aliscore/Download/index.en.html) - scripts tested with ALISCORE_v2.0
- RAxML - available from [here](https://github.com/stamatak/standard-RAxML) (required by process_genes.pl) - scripts tested with RAxML 7.5.4 
- ProteinModelSelection.pl - ships with RAxML (in usefulScripts directory) (required by process_genes.pl)


Walkthrough
-----------

The following section takes you step by step through the data preparation described in [Hahn et al. 2014](http://gbe.oxfordjournals.org/content/early/2014/04/13/gbe.evu078.short?rss=1 "Hahn et al. 2014 at GBE") - from the OrthoMCL output (groups.txt in data/ directory) to final concatenated alignment of genes with average bootstrap support >= 90:


	-bash-4.1$ select_clusters_v2.pl --groups groups.txt.gz --fasta goodProteins.fasta --critical critical.txt --exclusive > sort_clusters.log

optionally distribute the fasta files across a number of directories. THe below command would distribute all *.fasta files in the current directory (specified by ".") into sub-directories, with 50 *.fasta files per directory. Most of my downstream scripts serially process all files in the directory. So distributing files into separate directories now, could allow to parallelize the following steps.  

	-bash-4.1$ distribute.pl . 50

With the next command we will do the following to every *.fasta file in the directory:
-align sequences using clustalo
-perform alignment trimming with Aliscore and Alicut
-find the best fitting model of evolution using ProteinModelSelection.pl script and RAxML
-infer ML tree using the best fitting model using RAxML (100 rapid boostrap pseudoreplicates)

	-bash-4.1$ process_genes.pl . --nomask --bootstrap 100


to be continued..

