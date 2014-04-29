This repository contains useful script for the preparation of data for phylogenomic analyses.

The scripts in this repository have been used in the analyses presented in [Hahn et al. 2014](http://gbe.oxfordjournals.org/content/early/2014/04/13/gbe.evu078.short?rss=1 "Hahn et al. 2014 at GBE"). They were tailored specifically for this study, but might be useful for other users as well. Obviously, I won't guarantee that they will work for your data. If you find them useful, I would appreciate if you would cite our [paper](http://gbe.oxfordjournals.org/content/early/2014/04/13/gbe.evu078.short?rss=1 "Hahn et al. 2014 at GBE"). 

If you have questions/comments: christoph.hahn@nhm.uio.no

A step by step guide through some parts of the process described in [Hahn et al. 2014](http://gbe.oxfordjournals.org/content/early/2014/04/13/gbe.evu078.short?rss=1 "Hahn et al. 2014 at GBE") is given below.


PREREQUISITES:
--------------

A number of external scripts/programs are expected to be in your path:

sort_contigs.pl - script available from [here](http://www.genome.ou.edu/informatics.html) (required by sort_clusters_v2.pl if you want to extract the sequences into a fasta file - see below)

fasta2phylip.pl - script available from [here](http://www.cs.utexas.edu/~mswenson/misc/useful_scripts/) (required by process_genes.pl)

clustal omega - available from [here](http://www.clustal.org/omega/) (required by process_genes.pl) - scripts were tested with clustal omega 1.1 

ALICUT - available from [here](http://zfmk.de/web/Forschung/Abteilungen/AG_Wgele/Software/Utilities/index.en.html)(required by process_genes.pl) - scripts were tested with ALICUT_v2.3

ALICSCORE - availalbe from [here](http://zfmk.de/web/Forschung/Abteilungen/AG_Wgele/Software/Aliscore/Download/index.en.html) - scripts tested with ALISCORE_v2.0

RAxML - available from [here]() 

ProteinModelSelection.pl - ships with RAxML (in usefulScripts directory)


WALKTHROUGH:
------------
	-bash-4.1$ select_clusters_v2.pl --groups groups.txt.gz --fasta goodProteins.fasta --critical critical.txt --exclusive > sort_clusters.log

optionally distribute the fasta files across a number of directories:

	-bash-4.1$ distribute.pl . 50


to be continued..

