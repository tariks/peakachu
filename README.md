> **_NOTE:_**  Peakachu (version>=1.1.2) now supports both [.hic](https://github.com/aidenlab/juicer/wiki/Data) and [.cool](https://cooler.readthedocs.io/en/latest/datamodel.html) formats.

# Introduction

Accurately predicting chromatin loops from genome-wide interaction matrices such as Hi-C data is critical to deepening our understanding of proper gene regulation. Current approaches are mainly focused on searching for statistically enriched dots on a genome-wide map. However, given the availability of orthogonal data types such as ChIA-PET, HiChIP, Capture Hi-C, and high-throughput imaging, a supervised learning approach could facilitate the discovery of a comprehensive set of chromatin interactions. Here, we present Peakachu, a Random Forest classification framework that predicts chromatin loops from genome-wide contact maps. We compare Peakachu with current enrichment-based approaches, and find that Peakachu identifies a unique set of short-range interactions. We show that our models perform well in different platforms, across different sequencing depths, and across different species. We apply this framework to predict chromatin loops in 56 Hi-C datasets, and release the results at the 3D Genome Browser.

# Citation

Salameh, T.J., Wang, X., Song, F. et al. A supervised learning framework for chromatin loop detection in genome-wide contact maps. Nat Commun 11, 3428 (2020). https://doi.org/10.1038/s41467-020-17239-9

# Installation

Peakachu requires Python3 and several scientific packages to run. It is best to first set up the environment using conda and then install Peakachu from PyPI:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda create -n peakachu cooler scikit-learn numba joblib=1.1.0
conda activate peakachu
pip install -U peakachu hic-straw==0.0.6 
```

Peakachu should now be installed as a command-line tool within the new environment. Options for all peakachu commands and sub-commands can be accessed with the -h option. 


```bash
peakachu -h
```

    usage: peakachu [-h] {train,score_chromosome,score_genome,depth,pool} ...

    Unveil Hi-C Anchors and Peaks.

    positional arguments:
      {train,score_chromosome,score_genome,depth,pool}
        train               Train RandomForest model per chromosome
        score_chromosome    Calculate interaction probability per pixel for a chromosome
        score_genome        Calculate interaction probability per pixel for the whole genome
        depth               Calculate the total number of intra-chromosomal chromatin contacts and select the most appropriate pre-trained model
                            for you.
        pool                Print centroid loci from score_genome/score_chromosome output

    options:
      -h, --help            show this help message and exit


# Example: predicting loops in GM12878 Hi-C

The following example will download an example cooler file containing the GM12878 Hi-C data at the 10kb resolution, train a series of models using H3K27ac HiChIP interactions, and then predict loops using the trained models.

## Data preparation

Peakachu requires the contact map to be a .cool file or a .hic file and any training input to be a text file in bedpe format. Example training data can be found at the [training-sets](https://github.com/tariks/peakachu/tree/master/training-sets) subfolder. Cooler files may be found at the [4DN data portal](https://data.4dnucleome.org/).

```bash
wget -O Rao2014-GM12878-MboI-allreps-filtered.10kb.cool -L https://dl.dropboxusercontent.com/s/9d3hqyy0ou0rwsz/Rao2014-GM12878-MboI-allreps-filtered.10kb.cool?dl=0
```

## Train a model and predict loops

It is always a good idea to call the help function immediately before entering a command:

```bash
peakachu train -h
```

    usage: peakachu train [-h] [-r RESOLUTION] [-p PATH] [--balance] [-b BEDPE] [-w WIDTH] [--nproc NPROC] [-O OUTPUT]

    options:
      -h, --help            show this help message and exit
      -r RESOLUTION, --resolution RESOLUTION
                            Resolution in bp (default 10000)
      -p PATH, --path PATH  Path to a .cool URI string or a .hic file.
      --balance             Whether or not using the ICE/KR-balanced matrix.
      -b BEDPE, --bedpe BEDPE
                            Path to the bedpe file containing positive training set.
      -w WIDTH, --width WIDTH
                            Number of bins added to center of window. default width=5 corresponds to 11x11 windows
      --nproc NPROC         Number of worker processes that will be allocated for training. (default 4)
      -O OUTPUT, --output OUTPUT
                            Folder path to store trained models.

```bash
peakachu train -r 10000 -p Rao2014-GM12878-MboI-allreps-filtered.10kb.cool --balance -O models -b gm12878.mumbach.h3k27ac-hichip.hg19.bedpe
```

This will train 23 random forest models, each labeled by a chromosome. The model for every chromosome
was trained using interactions from all the other 22 chromosomes in the provided bedpe file. The purpose of this is to avoid Peakachu to predict loops from the same map it used for training, without overfitting. To use these models, you may either use the score_chromosome function to predict loops in only one chromosome, or the score_genome function to perform a genome-wide prediction.


```bash
peakachu score_chromosome -h
```

    usage: peakachu score_chromosome [-h] [-r RESOLUTION] [-p PATH] [--balance] [-C CHROM] [-m MODEL] [-l LOWER] [-u UPPER]
                                 [--minimum-prob MINIMUM_PROB] [-O OUTPUT]

    options:
      -h, --help            show this help message and exit
      -r RESOLUTION, --resolution RESOLUTION
                            Resolution in bp (default 10000)
      -p PATH, --path PATH  Path to a .cool URI string or a .hic file.
      --balance             Whether or not using the ICE/KR-balanced matrix.
      -C CHROM, --chrom CHROM
                            Chromosome label. Only contact data within the specified chromosome will be considered.
      -m MODEL, --model MODEL
                            Path to pickled model file.
      -l LOWER, --lower LOWER
                            Lower bound of distance between loci in bins (default 6).
      -u UPPER, --upper UPPER
                            Upper bound of distance between loci in bins (default 300).
      --minimum-prob MINIMUM_PROB
                            Only output pixels with probability score greater than this value (default 0.5)
      -O OUTPUT, --output OUTPUT
                            Output file name.

```bash
peakachu score_chromosome -r 10000 -p Rao2014-GM12878-MboI-allreps-filtered.10kb.cool --balance -O GM12878-chr2-scores.bedpe -C chr2 -m models/chr2.pkl 
peakachu pool -r 10000 -i GM12878-chr2-scores.bedpe -o GM12878-chr2-loops.bedpe -t .9
```

The pool function serves to select the most significant non-redundant results from per-pixel probabilities calculated by the score functions. It is recommended to try different probability thresholds to achieve the best sensitivity-specificity tradeoff. The output is a standard bedpe file with the 7th and the final column containing the predicted probability from the random forest model and the interaction frequency extracted from the contact matrix, respectively, to support further filtering. The results can be visualized in [juicebox](https://github.com/aidenlab/Juicebox) or [higlass](https://docs.higlass.io) by loading as 2D annotations. Here is an example screenshot of predicted GM12878 loops in juicer:
![Predicted loops from model trained on H3K27ac HiChIP interactions](https://github.com/tariks/peakachu/blob/master/example/gm12878-h3k27ac-loops.png)

# Using Peakachu as a standard loop caller

Models for predicting loops in Hi-C have been trained using CTCF ChIA-PET interactions, H3K27ac HiChIP interactions, and a high-confidence loop set (loops that can be detected by at least two orthogonal methods from CTCF ChIA-PET, Pol2 ChIA-PET, Hi-C, CTCF HiChIP, H3K27ac HiChIP, SMC1A HiChIP, H3K4me3 PLAC-Seq, and TrAC-Loop) as positive training samples, at a variety of read depths. Simply download the appropriate model file and directly run the score_genome/score_chromosome function if you want to detect chromatin loops on your own Hi-C or Micro-C maps.

If you are using Peakachu>=2.0, please select a model from the following table:

| Total intra reads | high-confidence (5kb)                                                                              | high-confidence (10kb)                                                                               | high-confidence (25kb)                                                                               |
| ----------------- | -------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------- |
| 2 billion         | [total 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.2billion.5kb.w6.pkl)   | [total 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.2billion.10kb.w6.pkl)   | [total 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.2billion.25kb.w5.pkl)   |
| 1.8 billion       | [90% 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.1.8billion.5kb.w6.pkl)   | [90% 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.1.8billion.10kb.w6.pkl)   | [90% 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.1.8billion.25kb.w5.pkl)   |
| 1.6 billion       | [80% 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.1.6billion.5kb.w6.pkl)   | [80% 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.1.6billion.10kb.w6.pkl)   | [80% 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.1.6billion.25kb.w5.pkl)   |
| 1.4 billion       | [70% 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.1.4billion.5kb.w6.pkl)   | [70% 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.1.4billion.10kb.w6.pkl)   | [70% 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.1.4billion.25kb.w5.pkl)   |
| 1.2 billion       | [60% 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.1.2billion.5kb.w6.pkl)   | [60% 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.1.2billion.10kb.w6.pkl)   | [60% 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.1.2billion.25kb.w5.pkl)   |
| 1 billion         | [50% 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.1billion.5kb.w6.pkl)     | [50% 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.1billion.10kb.w6.pkl)     | [50% 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.1billion.25kb.w5.pkl)     |
| 900 million       | [45% 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.900million.5kb.w6.pkl)   | [45% 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.900million.10kb.w6.pkl)   | [45% 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.900million.25kb.w5.pkl)   |
| 850 million       | [42.5% 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.850million.5kb.w6.pkl) | [42.5% 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.850million.10kb.w6.pkl) | [42.5% 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.850million.25kb.w5.pkl) |
| 800 million       | [40% 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.800million.5kb.w6.pkl)   | [40% 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.800million.10kb.w6.pkl)   | [40% 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.800million.25kb.w5.pkl)   |
| 750 million       | [37.5% 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.750million.5kb.w6.pkl) | [37.5% 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.750million.10kb.w6.pkl) | [37.5% 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.750million.25kb.w5.pkl) |
| 700 million       | [35% 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.700million.5kb.w6.pkl)   | [35% 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.700million.10kb.w6.pkl)   | [35% 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.700million.25kb.w5.pkl)   |
| 650 million       | [32.5% 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.650million.5kb.w6.pkl) | [32.5% 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.650million.10kb.w6.pkl) | [32.5% 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.650million.25kb.w5.pkl) |
| 600 million       | [30% 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.600million.5kb.w6.pkl)   | [30% 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.600million.10kb.w6.pkl)   | [30% 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.600million.25kb.w5.pkl)   |
| 550 million       | [27.5% 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.550million.5kb.w6.pkl) | [27.5% 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.550million.10kb.w6.pkl) | [27.5% 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.550million.25kb.w5.pkl) |
| 500 million       | [25% 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.500million.5kb.w6.pkl)   | [25% 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.500million.10kb.w6.pkl)   | [25% 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.500million.25kb.w5.pkl)   |
| 450 million       | [22.5% 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.450million.5kb.w6.pkl) | [22.5% 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.450million.10kb.w6.pkl) | [22.5% 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.450million.25kb.w5.pkl) |
| 400 million       | [20% 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.400million.5kb.w6.pkl)   | [20% 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.400million.10kb.w6.pkl)   | [20% 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.400million.25kb.w5.pkl)   |
| 350 million       | [17.5% 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.350million.5kb.w6.pkl) | [17.5% 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.350million.10kb.w6.pkl) | [17.5% 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.350million.25kb.w5.pkl) |
| 300 million       | [15% 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.300million.5kb.w6.pkl)   | [15% 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.300million.10kb.w6.pkl)   | [15% 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.300million.25kb.w5.pkl)   |
| 250 million       | [12.5% 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.250million.5kb.w6.pkl) | [12.5% 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.250million.10kb.w6.pkl) | [12.5% 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.250million.25kb.w5.pkl) |
| 200 million       | [10% 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.200million.5kb.w6.pkl)   | [10% 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.200million.10kb.w6.pkl)   | [10% 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.200million.25kb.w5.pkl)   |
| 150 million       | [7.5% 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.150million.5kb.w6.pkl)  | [7.5% 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.150million.10kb.w6.pkl)  | [7.5% 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.150million.25kb.w5.pkl)  |
| 100 million       | [5% 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.100million.5kb.w6.pkl)    | [5% 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.100million.10kb.w6.pkl)    | [5% 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.100million.25kb.w5.pkl)    |
| 50 million        | [2.5% 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.50million.5kb.w6.pkl)   | [2.5% 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.50million.10kb.w6.pkl)   | [2.5% 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.50million.25kb.w5.pkl)   |
| 30 million        | [1.5% 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.30million.5kb.w6.pkl)   | [1.5% 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.30million.10kb.w6.pkl)   | [1.5% 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.30million.25kb.w5.pkl)   |
| 10 million        | [0.5% 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.10million.5kb.w6.pkl)   | [0.5% 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.10million.10kb.w6.pkl)   | [0.5% 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.10million.25kb.w5.pkl)   |
| 5 million         | [0.25% 5kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.5million.5kb.w6.pkl)   | [0.25% 10kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.5million.10kb.w6.pkl)   | [0.25% 25kb](http://3dgenome.fsm.northwestern.edu/peakachu/high-confidence.5million.25kb.w5.pkl)   |

Instead, if you are using an older Peakachu version (<2.0), please select a model
from this table:

| Total intra reads | high-confidence (5kb)                                                                                    | high-confidence (10kb)                                                                                     | high-confidence (25kb)                                                                                     | CTCF Models (10kb)                                                                      | H3K27ac Model (10kb)                                                                          |
|-------------------|----------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------|
| 2 billion         | [total 5kb](https://dl.dropboxusercontent.com/s/4iq9lw3rta36oa7/high-confidence.2billion.5kb.pkl?dl=0)   | [total 10kb](https://dl.dropboxusercontent.com/s/3b8txqo5jj8cwga/high-confidence.2billion.10kb.pkl?dl=0)   | [total 25kb](https://dl.dropboxusercontent.com/s/l5wg6nmzd1t5nv5/high-confidence.2billion.25kb.pkl?dl=0)   | [CTCF total](https://dl.dropboxusercontent.com/s/enyg2m7ebj8mxsv/down100.ctcf.pkl?dl=0) | [H3K27ac total](https://dl.dropboxusercontent.com/s/yasl5hu0v510k2v/down100.h3k27ac.pkl?dl=0) |
| 1.8 billion       | [90% 5kb](https://dl.dropboxusercontent.com/s/vfp1pl424kwxbpp/high-confidence.1.8billion.5kb.pkl?dl=0)   | [90% 10kb](https://dl.dropboxusercontent.com/s/kuzzfqcwzso7lu0/high-confidence.1.8billion.10kb.pkl?dl=0)   | [90% 25kb](https://dl.dropboxusercontent.com/s/16tmjrq329uhzwt/high-confidence.1.8billion.25kb.pkl?dl=0)   | [CTCF 90%](https://dl.dropboxusercontent.com/s/g12hy9f28igh0ng/down90.ctcf.pkl?dl=0)    | [H3K27ac 90%](https://dl.dropboxusercontent.com/s/kdbv52eeilkzqfr/down90.h3k27ac.pkl?dl=0)    |
| 1.6 billion       | [80% 5kb](https://dl.dropboxusercontent.com/s/zgwbpt1rhhkirss/high-confidence.1.6billion.5kb.pkl?dl=0)   | [80% 10kb](https://dl.dropboxusercontent.com/s/wtkmw7zzvxsv7l2/high-confidence.1.6billion.10kb.pkl?dl=0)   | [80% 10kb](https://dl.dropboxusercontent.com/s/7tern39gqeph6tr/high-confidence.1.6billion.25kb.pkl?dl=0)   | [CTCF 80%](https://dl.dropboxusercontent.com/s/n2m4jxxojh0u5ay/down80.ctcf.pkl?dl=0)    | [H3K27ac 80%](https://dl.dropboxusercontent.com/s/45ekayzigeyuown/down80.h3k27ac.pkl?dl=0)    |
| 1.4 billion       | [70% 5kb](https://dl.dropboxusercontent.com/s/1iyxfbj3dn54l4y/high-confidence.1.4billion.5kb.pkl?dl=0)   | [70% 10kb](https://dl.dropboxusercontent.com/s/dqc8dcmkfexjdp8/high-confidence.1.4billion.10kb.pkl?dl=0)   | [70% 25kb](https://dl.dropboxusercontent.com/s/l55vokib8apx1bb/high-confidence.1.4billion.25kb.pkl?dl=0)   | [CTCF 70%](https://dl.dropboxusercontent.com/s/h9vm8z0uysti8xm/down70.ctcf.pkl?dl=0)    | [H3K27ac 70%](https://dl.dropboxusercontent.com/s/mrhe0uayv402vfk/down70.h3k27ac.pkl?dl=0)    |
| 1.2 billion       | [60% 5kb](https://dl.dropboxusercontent.com/s/5enu57up2qf1p2f/high-confidence.1.2billion.5kb.pkl?dl=0)   | [60% 10kb](https://dl.dropboxusercontent.com/s/684nnbu707p5pa9/high-confidence.1.2billion.10kb.pkl?dl=0)   | [60% 25kb](https://dl.dropboxusercontent.com/s/qgggtalcgf1islv/high-confidence.1.2billion.25kb.pkl?dl=0)   | [CTCF 60%](https://dl.dropboxusercontent.com/s/cfkfem4w8dhhgwm/down60.ctcf.pkl?dl=0)    | [H3K27ac 60%](https://dl.dropboxusercontent.com/s/0f9xv6ljjlcwnsv/down60.h3k27ac.pkl?dl=0)    |
| 1 billion         | [50% 5kb](https://dl.dropboxusercontent.com/s/gf10pbwd0fh5uvn/high-confidence.1billion.5kb.pkl?dl=0)     | [50% 10kb](https://dl.dropboxusercontent.com/s/crkmej6fvzwls82/high-confidence.1billion.10kb.pkl?dl=0)     | [50% 25kb](https://dl.dropboxusercontent.com/s/m5oyq2n5pajkic2/high-confidence.1billion.25kb.pkl?dl=0)     | [CTCF 50%](https://dl.dropboxusercontent.com/s/c0b6axxb16p2nd7/down50.ctcf.pkl?dl=0)    | [H3K27ac 50%](https://dl.dropboxusercontent.com/s/3w4befpvu7c7cqe/down50.h3k27ac.pkl?dl=0)    |
| 900 million       | [45% 5kb](https://dl.dropboxusercontent.com/s/z3tbnl8mvvaxu6a/high-confidence.900million.5kb.pkl?dl=0)   | [45% 10kb](https://dl.dropboxusercontent.com/s/ezdu69c0j6h0dql/high-confidence.900million.10kb.pkl?dl=0)   | [45% 25kb](https://dl.dropboxusercontent.com/s/3t3jsbsfeshkpmq/high-confidence.900million.25kb.pkl?dl=0)   |                                                                                         |                                                                                               |
| 850 million       | [42.5% 5kb](https://dl.dropboxusercontent.com/s/au1f2s15jgdqp5a/high-confidence.850million.5kb.pkl?dl=0) | [42.5% 10kb](https://dl.dropboxusercontent.com/s/4ymigubbmeaf341/high-confidence.850million.10kb.pkl?dl=0) | [42.5% 25kb](https://dl.dropboxusercontent.com/s/23jnf280qemda2f/high-confidence.850million.25kb.pkl?dl=0) |                                                                                         |                                                                                               |
| 800 million       | [40% 5kb](https://dl.dropboxusercontent.com/s/r7eciv6op3qg0ij/high-confidence.800million.5kb.pkl?dl=0)   | [40% 10kb](https://dl.dropboxusercontent.com/s/zr2tc8x7n615j2i/high-confidence.800million.10kb.pkl?dl=0)   | [40% 25kb](https://dl.dropboxusercontent.com/s/ape2vlf6rmomp3g/high-confidence.800million.25kb.pkl?dl=0)   | [CTCF 40%](https://dl.dropboxusercontent.com/s/8lvcdjenyoc8ggy/down40.ctcf.pkl?dl=0)    | [H3K27ac 40%](https://dl.dropboxusercontent.com/s/xwlk864nkoafzsy/down40.h3k27ac.pkl?dl=0)    |
| 750 million       | [37.5% 5kb](https://dl.dropboxusercontent.com/s/m7jah602d34xqqy/high-confidence.750million.5kb.pkl?dl=0) | [37.5% 10kb](https://dl.dropboxusercontent.com/s/cuw30u15a6vmla8/high-confidence.750million.10kb.pkl?dl=0) | [37.5% 25kb](https://dl.dropboxusercontent.com/s/f5gt3usivr8stqr/high-confidence.750million.25kb.pkl?dl=0) |                                                                                         |                                                                                               |
| 700 million       | [35% 5kb](https://dl.dropboxusercontent.com/s/oy615sddva5r45e/high-confidence.700million.5kb.pkl?dl=0)   | [35% 10kb](https://dl.dropboxusercontent.com/s/tlhjkx9zc3mg040/high-confidence.700million.10kb.pkl?dl=0)   | [35% 25kb](https://dl.dropboxusercontent.com/s/zbcvxszjkp7rye9/high-confidence.700million.25kb.pkl?dl=0)   |                                                                                         |                                                                                               |
| 650 million       | [32.5% 5kb](https://dl.dropboxusercontent.com/s/8oua8pwxrbemjp9/high-confidence.650million.5kb.pkl?dl=0) | [32.5% 10kb](https://dl.dropboxusercontent.com/s/0l5ilp0xodb8odh/high-confidence.650million.10kb.pkl?dl=0) | [32.5% 25kb](https://dl.dropboxusercontent.com/s/ofdkxau6sz1cshs/high-confidence.650million.25kb.pkl?dl=0) |                                                                                         |                                                                                               |
| 600 million       | [30% 5kb](https://dl.dropboxusercontent.com/s/p0nxkg8tstnzkoo/high-confidence.600million.5kb.pkl?dl=0)   | [30% 10kb](https://dl.dropboxusercontent.com/s/z3f0vx6t57a9cjg/high-confidence.600million.10kb.pkl?dl=0)   | [30% 25kb](https://dl.dropboxusercontent.com/s/gfolpr069uf3rnx/high-confidence.600million.25kb.pkl?dl=0)   | [CTCF 30%](https://dl.dropboxusercontent.com/s/f1383jpzj3addi4/down30.ctcf.pkl?dl=0)    | [H3K27ac 30%](https://dl.dropboxusercontent.com/s/dyvtyqvu3wpq3a5/down30.h3k27ac.pkl?dl=0)    |
| 550 million       | [27.5% 5kb](https://dl.dropboxusercontent.com/s/fqw6kbihm0jx4s4/high-confidence.550million.5kb.pkl?dl=0) | [27.5% 10kb](https://dl.dropboxusercontent.com/s/ymhq3mzw694z0xd/high-confidence.550million.10kb.pkl?dl=0) | [27.5% 25kb](https://dl.dropboxusercontent.com/s/sc6wdtkq6e1uav9/high-confidence.550million.25kb.pkl?dl=0) |                                                                                         |                                                                                               |
| 500 million       | [25% 5kb](https://dl.dropboxusercontent.com/s/4q1c3l1hghvladw/high-confidence.500million.5kb.pkl?dl=0)   | [25% 10kb](https://dl.dropboxusercontent.com/s/isy8jc90hlcrw3g/high-confidence.500million.10kb.pkl?dl=0)   | [25% 25kb](https://dl.dropboxusercontent.com/s/h9r980qx5cj4hbf/high-confidence.500million.25kb.pkl?dl=0)   |                                                                                         |                                                                                               |
| 450 million       | [22.5% 5kb](https://dl.dropboxusercontent.com/s/46701axrbjuuag2/high-confidence.450million.5kb.pkl?dl=0) | [22.5% 10kb](https://dl.dropboxusercontent.com/s/435z2ev5w6ldibz/high-confidence.450million.10kb.pkl?dl=0) | [22.5% 25kb](https://dl.dropboxusercontent.com/s/x8b69lusyprxmzn/high-confidence.450million.25kb.pkl?dl=0) |                                                                                         |                                                                                               |
| 400 million       | [20% 5kb](https://dl.dropboxusercontent.com/s/nrraex21xhmr2sk/high-confidence.400million.5kb.pkl?dl=0)   | [20% 10kb](https://dl.dropboxusercontent.com/s/jix51tru4emz1wl/high-confidence.400million.10kb.pkl?dl=0)   | [20% 25kb](https://dl.dropboxusercontent.com/s/mik9oe4b3ftj0e8/high-confidence.400million.25kb.pkl?dl=0)   | [CTCF 20%](https://dl.dropboxusercontent.com/s/a5nwa1xlg22ud24/down20.ctcf.pkl?dl=0)    | [H3K27ac 20%](https://dl.dropboxusercontent.com/s/qjm84cpw3uzlidp/down20.h3k27ac.pkl?dl=0)    |
| 350 million       | [17.5% 5kb](https://dl.dropboxusercontent.com/s/jpuvzlfkavx6aar/high-confidence.350million.5kb.pkl?dl=0) | [17.5% 10kb](https://dl.dropboxusercontent.com/s/qprrsdz0xsmfwzx/high-confidence.350million.10kb.pkl?dl=0) | [17.5% 25kb](https://dl.dropboxusercontent.com/s/hjefbyc81uevouh/high-confidence.350million.25kb.pkl?dl=0) |                                                                                         |                                                                                               |
| 300 million       | [15% 5kb](https://dl.dropboxusercontent.com/s/u60oaa9wqz2d5xj/high-confidence.300million.5kb.pkl?dl=0)   | [15% 10kb](https://dl.dropboxusercontent.com/s/i16s80swjtk0c7s/high-confidence.300million.10kb.pkl?dl=0)   | [15% 25kb](https://dl.dropboxusercontent.com/s/j1pj6e22w60x8f0/high-confidence.300million.25kb.pkl?dl=0)   |                                                                                         |                                                                                               |
| 250 million       | [12.5% 5kb](https://dl.dropboxusercontent.com/s/5if19vpjgfobght/high-confidence.250million.5kb.pkl?dl=0) | [12.5% 10kb](https://dl.dropboxusercontent.com/s/1winvhp5hcf7zer/high-confidence.250million.10kb.pkl?dl=0) | [12.5% 25kb](https://dl.dropboxusercontent.com/s/8bokmq6rpswi5yw/high-confidence.250million.25kb.pkl?dl=0) |                                                                                         |                                                                                               |
| 200 million       | [10% 5kb](https://dl.dropboxusercontent.com/s/i2axj1ij4vbhcha/high-confidence.200million.5kb.pkl?dl=0)   | [10% 10kb](https://dl.dropboxusercontent.com/s/joww6ej582rxz4f/high-confidence.200million.10kb.pkl?dl=0)   | [10% 25kb](https://dl.dropboxusercontent.com/s/in6v159h9w9i97p/high-confidence.200million.25kb.pkl?dl=0)   | [CTCF 10%](https://dl.dropboxusercontent.com/s/cqi0ws8een9ad4t/down10.ctcf.pkl?dl=0)    | [H3K27ac 10%](https://dl.dropboxusercontent.com/s/q8mlwn4mz6rnumr/down10.h3k27ac.pkl?dl=0)    |
| 150 million       | [7.5% 5kb](https://dl.dropboxusercontent.com/s/80qxv0jnrhdp1ud/high-confidence.150million.5kb.pkl?dl=0)  | [7.5% 10kb](https://dl.dropboxusercontent.com/s/ff8i1hjp6949rqy/high-confidence.150million.10kb.pkl?dl=0)  | [7.5% 25kb](https://dl.dropboxusercontent.com/s/k4yi68pj2wpg1d5/high-confidence.150million.25kb.pkl?dl=0)  |                                                                                         |                                                                                               |
| 100 million       | [5% 5kb](https://dl.dropboxusercontent.com/s/mwcdn9j4rjrb5u2/high-confidence.100million.5kb.pkl?dl=0)    | [5% 10kb](https://dl.dropboxusercontent.com/s/qw4ay9u6khh36bn/high-confidence.100million.10kb.pkl?dl=0)    | [5% 25kb](https://dl.dropboxusercontent.com/s/a7lf3p7pzfrts1c/high-confidence.100million.25kb.pkl?dl=0)    |                                                                                         |                                                                                               |
| 50 million        | [2.5% 5kb](https://dl.dropboxusercontent.com/s/5nnxf9drw39acu5/high-confidence.50million.5kb.pkl?dl=0)   | [2.5% 10kb](https://dl.dropboxusercontent.com/s/aeeu6r4wku5fe9s/high-confidence.50million.10kb.pkl?dl=0)   | [2.5% 25kb](https://dl.dropboxusercontent.com/s/rftwrdm1lq8inkd/high-confidence.50million.25kb.pkl?dl=0)   |                                                                                         |                                                                                               |
| 30 million        | [1.5% 5kb](https://dl.dropboxusercontent.com/s/3r653c7rioxf9ct/high-confidence.30million.5kb.pkl?dl=0)   | [1.5% 10kb](https://dl.dropboxusercontent.com/s/4mugruqjnx6tdul/high-confidence.30million.10kb.pkl?dl=0)   | [1.5% 25kb](https://dl.dropboxusercontent.com/s/g6nq0i0e44w33vi/high-confidence.30million.25kb.pkl?dl=0)   | [CTCF 1.5%](https://dl.dropboxusercontent.com/s/5gxeervadlga1b3/down1.ctcf.pkl?dl=0)    | [H3K27ac 1.5%](https://dl.dropboxusercontent.com/s/uh98lt1rbyauhgn/down1.h3k27ac.pkl?dl=0)    |
| 10 million        |                                                                                                          |                                                                                                            | [0.5% 25kb](https://dl.dropboxusercontent.com/s/gsaxkgz0oh4ahgf/high-confidence.10million.25kb.pkl?dl=0)   |                                                                                         |                                                                                               |
| 5 million         |                                                                                                          |                                                                                                            | [0.25% 25kb](https://dl.dropboxusercontent.com/s/10fbe85lpfabupw/high-confidence.5million.25kb.pkl?dl=0)   |                                                                                         |                                                                                               |

To make it clear, let's download another Hi-C dataset from 4DN: https://data.4dnucleome.org/files-processed/4DNFI5IHU27G/@@download/4DNFI5IHU27G.mcool. Peakachu provides a handy function `peakachu depth` to extract the total number of intra-chromosomal pairs in your data and help you select the most appropriate pre-trained model:


```bash
peakachu depth -p 4DNFI5IHU27G.mcool::resolutions/1000000
```

The output of above command will be:

    num of intra reads in your data: 592991890
    num of intra reads in a human with matched sequencing coverage: 582003409
    suggested model: 600 million

Therefore, we recommend using the 30% models (trained with ~600 million intra reads)
to predict loops on this data.

```bash
peakachu score_genome -r 10000 --balance -p 4DNFI5IHU27G.mcool::resolutions/10000 -O 4DNFI5IHU27G-peakachu-10kb-scores.bedpe -m high-confidence.600million.10kb.w6.pkl
peakachu pool -r 10000 -i 4DNFI5IHU27G-peakachu-10kb-scores.bedpe -o 4DNFI5IHU27G-peakachu-10kb-loops.0.95.bedpe -t 0.95
```

# Not just Hi-C
Peakachu has been tested on Hi-C, Micrco-C, and DNA SPRITE contact maps with good results. For training sets, ChIA-PET, HiChIP, and PLAC-Seq have been tested. The purpose of this software is ultimately to facilitate the interpretation of results from multiple types of experiments, and the user is encouraged to apply Peakachu's training framework to newer approaches as they become available.

# Release Notes
### Version 2.0 (09/06/2022)
1. Re-trained the models using the latest scikit-learn v1.1.2
2. Used the distance-normalized signals instead of original contact signals
3. Added a 2D Gaussian filter followed by min-max scaling to pre-process each training image
4. Optimized the computation efficiency using numba and matrix operations.