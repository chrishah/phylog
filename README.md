Introduction
------------

This repository contains useful script for the preparation of data for phylogenomic analyses.

The scripts in this repository were used in the analyses presented in [Hahn et al. 2014](http://gbe.oxfordjournals.org/content/early/2014/04/13/gbe.evu078.short?rss=1 "Hahn et al. 2014 at GBE"). They were tailored specifically for this study, but might be useful for other users as well. Obviously, I won't guarantee that they will work for your data. If you find them useful, I would appreciate if you would cite our [paper](http://gbe.oxfordjournals.org/content/early/2014/04/13/gbe.evu078.short?rss=1 "Hahn et al. 2014 at GBE"). 

If you have questions/comments don't hesitate to contact me at: c.hahn@hull.ac.uk

A step by step guide through some parts of the process described in [Hahn et al. 2014](http://gbe.oxfordjournals.org/content/early/2014/04/13/gbe.evu078.short?rss=1 "Hahn et al. 2014 at GBE") is given below. For references to the external programs used by my scripts please refer to the methods section of our paper.


Prerequisites
-------------

A number of external scripts/programs are expected to be in your path:

- sort_contigs.pl - script available from [here](http://www.genome.ou.edu/informatics.html) (required by sort_clusters_v2.pl if you want to extract the sequences into a fasta file - see below)
- fasta2phylip.pl - script available from [here](https://github.com/chinchliff/physcripts/blob/master/fasta2phylip.pl) (required by process_genes.pl)
- clustal omega - available from [here](http://www.clustal.org/omega/) (required by process_genes.pl) - scripts were tested with clustal omega 1.1 
- ALICUT - available from [here](https://www.zfmk.de/en/research/research-centres-and-groups/utilities) (required by process_genes.pl) - scripts were tested with ALICUT_v2.3
- ALICSCORE - availalbe from [here](https://www.zfmk.de/en/research/research-centres-and-groups/aliscore) (required by process_genes.pl) - scripts tested with Aliscore.02.2.pl
- RAxML - available from [here](https://github.com/stamatak/standard-RAxML) (required by process_genes.pl) - scripts tested with RAxML 7.5.4 
- ProteinModelSelection.pl - ships with RAxML (in usefulScripts directory) (required by process_genes.pl)


Walkthrough
-----------

The following section takes you step by step through the data preparation described in [Hahn et al. 2014](http://gbe.oxfordjournals.org/content/early/2014/04/13/gbe.evu078.short?rss=1 "Hahn et al. 2014 at GBE") - from the OrthoMCL output (groups.txt in data/ directory) to final concatenated alignment of genes with average bootstrap support >= 90:

__We start__ by fitering the orthology clustering result obtained by OrthoMCL. The below command will apply the following filters:
+ minimal number of taxa contained in the cluster, default: 8
+ maximum median number of sequences per taxon, default: 2
+ maximum mean number of sequences per taxon, default: 5
+ cluster contains all taxa listed in critical.txt
+ only the taxa listed in critical.txt will be retained

A fasta file containing the protein sequences for each of the retained cluster will be written to `cluster_select_output`.
See a list of all options by simply running `select_clusters_v2.pl` without options.
```bash
-bash-4.1$ select_clusters_v2.pl --groups groups.txt.gz --fasta goodProteins.fasta --critical critical.txt --exclusive > sort_clusters.log
```
__Next__ we will infer ML trees for each individual alignment. _Optionally_ distribute the resulting fasta files across a number of directories. The below command would distribute all *.fasta files in the current directory (specified by ".") into sub-directories, with 50 *.fasta files per directory. Most of my downstream scripts serially process all files in a directory. To distribute files into separate directories may allow parallel processing of the directory content, if the required computational resources are available.
```bash
-bash-4.1$ cd cluster_select_output
-bash-4.1$ mkdir 1-initial-trees
-bash-4.1$ cd 1-initial-trees
-bash-4.1$ mv ../*.fasta .
-bash-4.1$ distribute.pl . 50
```
With the next command we will do the following to every *.fasta file in the directory `10000` (:
+ align sequences using clustalo
+ perform alignment trimming with Aliscore and Alicut
+ find the best fitting model of evolution using ProteinModelSelection.pl script and RAxML
+ infer ML tree using the best fitting model using RAxML (100 rapid boostrap pseudoreplicates)
Again, a list of all flags for the script is displayed when running `process_genes.pl` without options. 
```bash
-bash-4.1$ cd 10000
-bash-4.1$ process_genes.pl . --nomask --bootstrap 100
```

to be continued..

