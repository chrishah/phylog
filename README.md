Introduction
------------

This repository contains useful script for the preparation of data for phylogenomic analyses.

The scripts in this repository were used in the analyses presented in [Hahn et al. 2014](http://gbe.oxfordjournals.org/content/early/2014/04/13/gbe.evu078.short?rss=1 "Hahn et al. 2014 at GBE"). They were tailored specifically for this study, but might be useful for other users as well. Obviously, I won't guarantee that they will work for your data. If you find them useful, I would appreciate if you would cite our [paper](http://gbe.oxfordjournals.org/content/early/2014/04/13/gbe.evu078.short?rss=1 "Hahn et al. 2014 at GBE"). 

If you have questions/comments don't hesitate to contact me at: c.hahn@hull.ac.uk

A step by step guide through some parts of the process described in [Hahn et al. 2014](http://gbe.oxfordjournals.org/content/early/2014/04/13/gbe.evu078.short?rss=1 "Hahn et al. 2014 at GBE") is given below. For references to the external programs used by my scripts please refer to the methods section of our paper.


Prerequisites
-------------

A number of external scripts/programs are expected to be in your path (I have collected some of the external scripts [here](https://github.com/chrishah/phylog/tree/master/scripts-external):

- sort_contigs.pl - script available from [here](http://www.genome.ou.edu/informatics.html) (required by sort_clusters_v2.pl if you want to extract the sequences into a fasta file - see below)
- select_contigs.pl - script available from [here](http://www.genome.ou.edu/informatics.html)
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

A fasta file containing the protein sequences for each of the retained cluster will be written to the directory `cluster_select_output`. The file `sort_clusters.log` will contain information on the decision making for all clusters and a summary of the number of clusters retained/dropped at the end.
See a list of all options by simply running `select_clusters_v2.pl` without options.
```bash
select_clusters_v2.pl --groups groups.txt.gz --fasta goodProteins.fasta --critical critical.txt --exclusive > sort_clusters.log
```
__Next__ we will infer ML trees for each individual alignment. _Optionally_ distribute the resulting fasta files across a number of directories. The below command would distribute all *.fasta files in the current directory (specified by ".") into sub-directories, with 50 *.fasta files per directory. Most of my downstream scripts serially process all files in a directory. To distribute files into separate directories may allow parallel processing of the directory content, if the required computational resources are available.
```bash
cd cluster_select_output
WORKING=$(pwd)
mkdir 1-initial-trees
cd 1-initial-trees
mv ../*.fasta .
distribute.pl . 50
```
With the next command we will do the following to every *.fasta file in a directory (e.g.: `10000` - repeat/extend for data in additional directories accordingly):
+ align sequences using clustalo
+ perform alignment trimming with Aliscore and Alicut
+ find the best fitting model of evolution using ProteinModelSelection.pl script and RAxML
+ infer ML tree using the best fitting model using RAxML (100 rapid boostrap pseudoreplicates)
Again, a list of all flags for the script is displayed when running `process_genes.pl` without options. 
```bash
process_genes.pl ./10000 --nomask --bootstrap 100
```
__Next__ we will bin the resulting trees into three categories:
+ trees containing no paralogs (each taxon is represented by only one sequence)
+ trees containing in-paralogs (i.e. all paralogs form monopyletic groups)
+ trees containing out-paralogs (i.e. paralogs do not form monophyletic groups)
```bash
cd ..
mkdir 2-bin-trees
cd 2-bin-trees
#create list of all trees
for dir in $(ls -1 | grep "^[0-9]\{5\}" | sed 's/\/$//g'); do ls -1 $WORKING/1-initial-trees/$dir/*_processed/RAxML_bipartitions.ALICUT_*; done |perl -ne 'chomp; @a=split/\//; @b=split("_",$a[-1]); $out = $b[2]."_".$b[3];print substr($out,0,-4). "\n";' > IDs.list
#create symbolic links to all trees
mkdir trees
for file in $(<IDs.list); do ln -sv $WORKING/1-initial-trees/100??/$file\_processed/RAxML_bipartitions.ALICUT_$file.aln trees/RAxML_bipartitions.ALICUT_$file.aln; done
#bin trees
tree_bin_v2.pl trees > tree_bin.log
```

to be continued..

