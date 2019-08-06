
# Introduction

## What is Peakachu

Peakachu is an acronym that standands for Unveil Hi-C Anchors and Peaks. It takes genome-wide contact data as input and returns coordinates of likely interactions such as chromatin loops. A machine learning framework based on sklearn is employed to generate random forest models trained on example interactions predicted by an arbitrary experiment. For example, loops can be predicted in Hi-C data using models trained with the Hi-C profiles of interactions detected via ChIA-PET. Although Hi-C is the namesake of Peakachu, it is designed to accept any genome-wide contact map including those from Micro-C and DNA SPRITE.

## Installation

Peakachu requires Python3 and several scientific packages to run. It is best to set up a conda environment then install from github. Copy and paste the commands below:


```bash
conda create -n 3dgenome python=3.6
conda install -c bioconda cooler
source activate 3dgenome
git clone https://github.com/tariks/peakachu
cd peakachu
python setup.py --install
```

Peakachu should now be installed as a command-line tool within the new environment. Options for all peakachu commands can be accessed with the -h option.


```bash
peakachu -h
```

    usage: peakachu [-h] {train,score_chromosome,score_genome,depth,pool} ...
    
    Train Random Forest with Hi-C data and training peaklist.
    
    positional arguments:
      {train,score_chromosome,score_genome,depth,pool}
        train               Train RandomForest model per chromosome
        score_chromosome    Calculate interaction probability per pixel for a
                            chromosome
        score_genome        Calculate interaction probability per pixel for the
                            whole genome
        depth               Print read count of a .cool file
        pool                Print centroid loci from score_x output
    
    optional arguments:
      -h, --help            show this help message and exit



```bash
peakachu train -h
```

    usage: peakachu train [-h] [-p PATH] [-O OUTPUT] [-w WIDTH] [-b BEDPE]
    
    optional arguments:
      -h, --help            show this help message and exit
      -p PATH, --path PATH  URI string pointing to a .cool or multi-res .cool
                            file. Append ::10000 to the filename of a multi-res
                            .cool to use the 10kb matrix.
      -O OUTPUT, --output OUTPUT
                            Folder path to store results.
      -w WIDTH, --width WIDTH
                            number of bins added to center of window. default
                            width=5 corresponds to 11x11 windows
      -b BEDPE, --bedpe BEDPE
                            Path to the bedpe file containing pairs of resolutions
                            corresponding to experimentally detected loops.


# Example: predicting loops in GM12878 Hi-C

GM12878 is a commonly studied cell-line based on lymphoblasts from an adult individual. The following example will download a cooler file from a public source, train a series of models using ChIA-PET or HiChIP data, then predict loops using the trained models.

## Data preparation

Peakachu requires the contact map to be a cooler file and any training input to be a text file in bedpe format. Example training data is included in the github repository, consisting of bedpe files prepared from supplementary tables in [Tang et al][1] and [Mumbach et al][2]. Cooler files may be found at the [4DN data portal][3]
[1]: https://www.cell.com/cell/fulltext/S0092-8674(15)01504-4
[2]: https://www.ncbi.nlm.nih.gov/pubmed/28945252
[3]: https://data.4dnucleome.org/


```bash
wget ftp://cooler.csail.mit.edu/coolers/hg19/Rao2014-GM12878-MboI-allreps-filtered.10kb.cool
```

## Predicting loops in GM12878

It is always a good idea to call the help function immediately before entering a command


```bash
peakachu train -h
```

    usage: peakachu train [-h] [-p PATH] [-O OUTPUT] [-w WIDTH] [-b BEDPE]
    
    optional arguments:
      -h, --help            show this help message and exit
      -p PATH, --path PATH  URI string pointing to a .cool or multi-res .cool
                            file. Append ::10000 to the filename of a multi-res
                            .cool to use the 10kb matrix.
      -O OUTPUT, --output OUTPUT
                            Folder path to store results.
      -w WIDTH, --width WIDTH
                            number of bins added to center of window. default
                            width=5 corresponds to 11x11 windows
      -b BEDPE, --bedpe BEDPE
                            Path to the bedpe file containing pairs of resolutions
                            corresponding to experimentally detected loops.



```bash
peakachu train -p Rao2014-GM12878-MboI-allreps-filtered.10kb.cool -O models -b hg19.mumbach.h3k27ac.hichip.bedpe 
```

This will train one 23 random forest models, each labeled by a chromosome. Every model was trained on all of the interactions from the bedpe files EXCEPT for the chromosome which it is labeled as. The purpose of this is to allow Peakachu to predict loops from the same map it used for training, without overfitting. To use these models, you may either use the score_chromosome function to predict loops in only one chromosome, or the score_genome function when using a trained model to predict loops in a new contact map.


```bash
peakachu score_chromosome -h
```

    usage: peakachu score_chromosome [-h] [-p PATH] [-O OUTPUT] [-w WIDTH]
                                     [-m MODEL] [-l LOWER] [-u UPPER]
    
    optional arguments:
      -h, --help            show this help message and exit
      -p PATH, --path PATH  URI string pointing to a .cool or multi-res .cool
                            file. Append ::10000 to the filename of a multi-res
                            .cool to use the 10kb matrix.
      -O OUTPUT, --output OUTPUT
                            Folder path to store results.
      -w WIDTH, --width WIDTH
                            number of bins added to center of window. default
                            width=5 corresponds to 11x11 windows
      -m MODEL, --model MODEL
                            Path to pickled model file.
      -l LOWER, --lower LOWER
                            Lower bound of distance between loci in bins (default
                            2).
      -u UPPER, --upper UPPER
                            Lower bound of distance between loci in bins (default
                            300).



```bash
for i in models/*pkl; do peakachu score_chromosome -O scores -m $i; done
for i in scores/*; do peakachu pool -i $i -t .9 > ${i}.loops.txt; done
```

The pool function serves to select the most significant non-redundant results from per-pixel probabilities calculated by the score functions. It is recommended to try different probability thresholds to achieve the best sensitivity-specificity tradeoff.
