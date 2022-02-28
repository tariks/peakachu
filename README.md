> **_NOTE:_**  Peakachu (version>=1.1.2) now supports both [.hic](https://github.com/aidenlab/juicer/wiki/Data) and [.cool](https://cooler.readthedocs.io/en/latest/datamodel.html) formats.

# Introduction
## What is Peakachu
Peakachu is an acronym that standands for Unveil Hi-C Anchors and Peaks. It takes genome-wide contact data as input and returns coordinates of likely interactions such as chromatin loops. A machine learning framework based on sklearn is employed to generate random forest models trained on example interactions predicted by an arbitrary experiment. For example, loops can be predicted in Hi-C data using models trained with the Hi-C profiles of interactions detected via ChIA-PET. Although Hi-C is the namesake of Peakachu, it is designed to accept any genome-wide contact map including those from Micro-C and DNA SPRITE.

## Citation
Salameh, T.J., Wang, X., Song, F. et al. A supervised learning framework for chromatin loop detection in genome-wide contact maps. Nat Commun 11, 3428 (2020). https://doi.org/10.1038/s41467-020-17239-9

## Installation
Peakachu requires Python3 and several scientific packages to run. It is best to set up a conda environment then install from github. Copy and paste the command snippets below:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda create -n 3dgenome python=3.6 scikit-learn=0.20.2 numpy scipy pandas h5py cooler
source activate 3dgenome
pip install hic-straw
git clone https://github.com/tariks/peakachu
cd peakachu
python setup.py install
```

Peakachu should now be installed as a command-line tool within the new environment. Options for all peakachu commands and sub-commands can be accessed with the -h option. 


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
        depth               Print total intra-chromosomal read count
        pool                Print centroid loci from score_genome/score_chromosome
                            output
    
    optional arguments:
      -h, --help            show this help message and exit


# Example: predicting loops in GM12878 Hi-C

The following example will download an example cooler file containing the GM12878 Hi-C data at the 10kb resolution, train a series of models using H3K27ac HiChIP interactions, then predict loops using the trained models.

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

    usage: peakachu train [-h] [-r RESOLUTION] [-p PATH] [--balance] [-w WIDTH]
                      [-O OUTPUT] [-b BEDPE] [--nproc NPROC]

    optional arguments:
      -h, --help            show this help message and exit
      -r RESOLUTION, --resolution RESOLUTION
                            Resolution in bp (default 10000)
      -p PATH, --path PATH  Path to a .cool URI string or a .hic file.
      --balance             Whether or not using the ICE/KR-balanced matrix.
      -w WIDTH, --width WIDTH
                            Number of bins added to center of window. default
                            width=5 corresponds to 11x11 windows
      -O OUTPUT, --output OUTPUT
                            Folder path to store trained models.
      -b BEDPE, --bedpe BEDPE
                            Path to the bedpe file containing positive training
                            set.
      --nproc NPROC         Number of worker processes that will be allocated for
                            training. (default 4)

```bash
peakachu train -r 10000 -p Rao2014-GM12878-MboI-allreps-filtered.10kb.cool --balance -O models -b gm12878.mumbach.h3k27ac-hichip.hg19.bedpe
```

This will train 23 random forest models, each labeled by a chromosome. The model for every chromosome
was trained using interactions from all the other 22 chromosomes in the provided bedpe file. The purpose of this is to avoid Peakachu to predict loops from the same map it used for training, without overfitting. To use these models, you may either use the score_chromosome function to predict loops in only one chromosome, or the score_genome function to perform a genome-wide prediction.


```bash
peakachu score_chromosome -h
```

    usage: peakachu score_chromosome [-h] [-r RESOLUTION] [-p PATH] [--balance]
                                 [-w WIDTH] [-O OUTPUT] [-C CHROM] [-m MODEL]
                                 [-l LOWER] [-u UPPER]
                                 [--minimum-prob MINIMUM_PROB]

    optional arguments:
      -h, --help            show this help message and exit
      -r RESOLUTION, --resolution RESOLUTION
                            Resolution in bp (default 10000)
      -p PATH, --path PATH  Path to a .cool URI string or a .hic file.
      --balance             Whether or not using the ICE/KR-balanced matrix.
      -w WIDTH, --width WIDTH
                            Number of bins added to center of window. default
                            width=5 corresponds to 11x11 windows
      -O OUTPUT, --output OUTPUT
                            Output file name.
      -C CHROM, --chrom CHROM
                            Chromosome label. Only contact data within the
                            specified chromosome will be considered.
      -m MODEL, --model MODEL
                            Path to pickled model file.
      -l LOWER, --lower LOWER
                            Lower bound of distance between loci in bins (default
                            6).
      -u UPPER, --upper UPPER
                            Upper bound of distance between loci in bins (default
                            300).
      --minimum-prob MINIMUM_PROB
                            Only output pixels with probability score greater than
                            this value (default 0.5)

```bash
peakachu score_chromosome -r 10000 -p Rao2014-GM12878-MboI-allreps-filtered.10kb.cool --balance -O GM12878-chr2-scores.bedpe -C chr2 -m models/chr2.pkl 
peakachu pool -r 10000 -i GM12878-chr2-scores.bedpe -o GM12878-chr2-loops.bedpe -t .9
```

The pool function serves to select the most significant non-redundant results from per-pixel probabilities calculated by the score functions. It is recommended to try different probability thresholds to achieve the best sensitivity-specificity tradeoff. The output is a standard bedpe file with the 7th and the final column containing the predicted probability from the random forest model and the interaction frequency extracted from the contact matrix, respectively, to support further filtering. The results can be visualized in [juicebox](https://github.com/aidenlab/Juicebox) or [higlass](https://docs.higlass.io) by loading as 2D annotations. Here is an example screenshot of predicted GM12878 loops in juicer:
![Predicted loops from model trained on H3K27ac HiChIP interactions](https://github.com/tariks/peakachu/blob/master/example/gm12878-h3k27ac-loops.png)

# Using Peakachu as a standard loop caller

Models for predicting loops in Hi-C have been trained using CTCF ChIA-PET interactions, H3K27ac HiChIP interactions, and a high-confidence loop set (loops that can be detected by at least two orthogonal methods from CTCF ChIA-PET, Pol2 ChIA-PET, Hi-C, CTCF HiChIP, H3K27ac HiChIP, SMC1A HiChIP, H3K4me3 PLAC-Seq, and TrAC-Loop) as positive training samples, at a variety of read depths. Simply download the appropriate model file and directly run the score_genome/score_chromosome function if you want to
detect chromatin loops on your own Hi-C or Micro-C maps.

| Total intra reads | high-confidence (5kb)                                                                                  | high-confidence (10kb)                                                                                   | high-confidence (25kb)                                                                                   | CTCF Models (10kb)                                                                      | H3K27ac Model (10kb)                                                                          |
|-------------------|--------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------|
| 2 billion         | [total 5kb](https://dl.dropboxusercontent.com/s/4iq9lw3rta36oa7/high-confidence.2billion.5kb.pkl?dl=0) | [total 10kb](https://dl.dropboxusercontent.com/s/3b8txqo5jj8cwga/high-confidence.2billion.10kb.pkl?dl=0) | [total 25kb](https://dl.dropboxusercontent.com/s/l5wg6nmzd1t5nv5/high-confidence.2billion.25kb.pkl?dl=0) | [CTCF total](https://dl.dropboxusercontent.com/s/enyg2m7ebj8mxsv/down100.ctcf.pkl?dl=0) | [H3K27ac total](https://dl.dropboxusercontent.com/s/yasl5hu0v510k2v/down100.h3k27ac.pkl?dl=0) |
| 1.8 billion       | [90% 5kb](https://dl.dropboxusercontent.com/s/vfp1pl424kwxbpp/high-confidence.1.8billion.5kb.pkl?dl=0) | [90% 10kb](https://dl.dropboxusercontent.com/s/kuzzfqcwzso7lu0/high-confidence.1.8billion.10kb.pkl?dl=0) | [90% 25kb](https://dl.dropboxusercontent.com/s/16tmjrq329uhzwt/high-confidence.1.8billion.25kb.pkl?dl=0) | [CTCF 90%](https://dl.dropboxusercontent.com/s/g12hy9f28igh0ng/down90.ctcf.pkl?dl=0)    | [H3K27ac 90%](https://dl.dropboxusercontent.com/s/kdbv52eeilkzqfr/down90.h3k27ac.pkl?dl=0)    |
| 1.6 billion       | [80% 5kb](https://dl.dropboxusercontent.com/s/zgwbpt1rhhkirss/high-confidence.1.6billion.5kb.pkl?dl=0) | [80% 10kb](https://dl.dropboxusercontent.com/s/wtkmw7zzvxsv7l2/high-confidence.1.6billion.10kb.pkl?dl=0) | [80% 10kb](https://dl.dropboxusercontent.com/s/7tern39gqeph6tr/high-confidence.1.6billion.25kb.pkl?dl=0) | [CTCF 80%](https://dl.dropboxusercontent.com/s/n2m4jxxojh0u5ay/down80.ctcf.pkl?dl=0)    | [H3K27ac 80%](https://dl.dropboxusercontent.com/s/45ekayzigeyuown/down80.h3k27ac.pkl?dl=0)    |
| 1.4 billion       | [70% 5kb](https://dl.dropboxusercontent.com/s/1iyxfbj3dn54l4y/high-confidence.1.4billion.5kb.pkl?dl=0) | [70% 10kb](https://dl.dropboxusercontent.com/s/dqc8dcmkfexjdp8/high-confidence.1.4billion.10kb.pkl?dl=0) | [70% 25kb](https://dl.dropboxusercontent.com/s/l55vokib8apx1bb/high-confidence.1.4billion.25kb.pkl?dl=0) | [CTCF 70%](https://dl.dropboxusercontent.com/s/h9vm8z0uysti8xm/down70.ctcf.pkl?dl=0)    | [H3K27ac 70%](https://dl.dropboxusercontent.com/s/mrhe0uayv402vfk/down70.h3k27ac.pkl?dl=0)    |
| 1.2 billion       | [60% 5kb](https://dl.dropboxusercontent.com/s/5enu57up2qf1p2f/high-confidence.1.2billion.5kb.pkl?dl=0) | [60% 10kb](https://dl.dropboxusercontent.com/s/684nnbu707p5pa9/high-confidence.1.2billion.10kb.pkl?dl=0) | [60% 25kb](https://dl.dropboxusercontent.com/s/qgggtalcgf1islv/high-confidence.1.2billion.25kb.pkl?dl=0) | [CTCF 60%](https://dl.dropboxusercontent.com/s/cfkfem4w8dhhgwm/down60.ctcf.pkl?dl=0)    | [H3K27ac 60%](https://dl.dropboxusercontent.com/s/0f9xv6ljjlcwnsv/down60.h3k27ac.pkl?dl=0)    |
| 1 billion         | [50% 5kb](https://dl.dropboxusercontent.com/s/gf10pbwd0fh5uvn/high-confidence.1billion.5kb.pkl?dl=0)   | [50% 10kb](https://dl.dropboxusercontent.com/s/crkmej6fvzwls82/high-confidence.1billion.10kb.pkl?dl=0)   | [50% 25kb](https://dl.dropboxusercontent.com/s/m5oyq2n5pajkic2/high-confidence.1billion.25kb.pkl?dl=0)   | [CTCF 50%](https://dl.dropboxusercontent.com/s/c0b6axxb16p2nd7/down50.ctcf.pkl?dl=0)    | [H3K27ac 50%](https://dl.dropboxusercontent.com/s/3w4befpvu7c7cqe/down50.h3k27ac.pkl?dl=0)    |
| 800 million       | [40% 5kb](https://dl.dropboxusercontent.com/s/r7eciv6op3qg0ij/high-confidence.800million.5kb.pkl?dl=0) | [40% 10kb](https://dl.dropboxusercontent.com/s/zr2tc8x7n615j2i/high-confidence.800million.10kb.pkl?dl=0) | [40% 25kb](https://dl.dropboxusercontent.com/s/ape2vlf6rmomp3g/high-confidence.800million.25kb.pkl?dl=0) | [CTCF 40%](https://dl.dropboxusercontent.com/s/8lvcdjenyoc8ggy/down40.ctcf.pkl?dl=0)    | [H3K27ac 40%](https://dl.dropboxusercontent.com/s/xwlk864nkoafzsy/down40.h3k27ac.pkl?dl=0)    |
| 600 million       | [30% 5kb](https://dl.dropboxusercontent.com/s/p0nxkg8tstnzkoo/high-confidence.600million.5kb.pkl?dl=0) | [30% 10kb](https://dl.dropboxusercontent.com/s/z3f0vx6t57a9cjg/high-confidence.600million.10kb.pkl?dl=0) | [30% 25kb](https://dl.dropboxusercontent.com/s/gfolpr069uf3rnx/high-confidence.600million.25kb.pkl?dl=0) | [CTCF 30%](https://dl.dropboxusercontent.com/s/f1383jpzj3addi4/down30.ctcf.pkl?dl=0)    | [H3K27ac 30%](https://dl.dropboxusercontent.com/s/dyvtyqvu3wpq3a5/down30.h3k27ac.pkl?dl=0)    |
| 400 million       | [20% 5kb](https://dl.dropboxusercontent.com/s/nrraex21xhmr2sk/high-confidence.400million.5kb.pkl?dl=0) | [20% 10kb](https://dl.dropboxusercontent.com/s/jix51tru4emz1wl/high-confidence.400million.10kb.pkl?dl=0) | [20% 25kb](https://dl.dropboxusercontent.com/s/mik9oe4b3ftj0e8/high-confidence.400million.25kb.pkl?dl=0) | [CTCF 20%](https://dl.dropboxusercontent.com/s/a5nwa1xlg22ud24/down20.ctcf.pkl?dl=0)    | [H3K27ac 20%](https://dl.dropboxusercontent.com/s/qjm84cpw3uzlidp/down20.h3k27ac.pkl?dl=0)    |
| 200 million       | [10% 5kb](https://dl.dropboxusercontent.com/s/i2axj1ij4vbhcha/high-confidence.200million.5kb.pkl?dl=0) | [10% 10kb](https://dl.dropboxusercontent.com/s/joww6ej582rxz4f/high-confidence.200million.10kb.pkl?dl=0) | [10% 25kb](https://dl.dropboxusercontent.com/s/in6v159h9w9i97p/high-confidence.200million.25kb.pkl?dl=0) | [CTCF 10%](https://dl.dropboxusercontent.com/s/cqi0ws8een9ad4t/down10.ctcf.pkl?dl=0)    | [H3K27ac 10%](https://dl.dropboxusercontent.com/s/q8mlwn4mz6rnumr/down10.h3k27ac.pkl?dl=0)    |
| 100 million       | [5% 5kb](https://dl.dropboxusercontent.com/s/mwcdn9j4rjrb5u2/high-confidence.100million.5kb.pkl?dl=0)  | [5% 10kb](https://dl.dropboxusercontent.com/s/qw4ay9u6khh36bn/high-confidence.100million.10kb.pkl?dl=0)  | [5% 25kb](https://dl.dropboxusercontent.com/s/a7lf3p7pzfrts1c/high-confidence.100million.25kb.pkl?dl=0)  |                                                                                         |                                                                                               |
| 50 million        | [2.5% 5kb](https://dl.dropboxusercontent.com/s/5nnxf9drw39acu5/high-confidence.50million.5kb.pkl?dl=0) | [2.5% 10kb](https://dl.dropboxusercontent.com/s/aeeu6r4wku5fe9s/high-confidence.50million.10kb.pkl?dl=0) | [2.5% 25kb](https://dl.dropboxusercontent.com/s/rftwrdm1lq8inkd/high-confidence.50million.25kb.pkl?dl=0) |                                                                                         |                                                                                               |
| 30 million        | [1.5% 5kb](https://dl.dropboxusercontent.com/s/3r653c7rioxf9ct/high-confidence.30million.5kb.pkl?dl=0) | [1.5% 10kb](https://dl.dropboxusercontent.com/s/4mugruqjnx6tdul/high-confidence.30million.10kb.pkl?dl=0) | [1.5% 25kb](https://dl.dropboxusercontent.com/s/g6nq0i0e44w33vi/high-confidence.30million.25kb.pkl?dl=0) | [CTCF 1.5%](https://dl.dropboxusercontent.com/s/5gxeervadlga1b3/down1.ctcf.pkl?dl=0)    | [H3K27ac 1.5%](https://dl.dropboxusercontent.com/s/uh98lt1rbyauhgn/down1.h3k27ac.pkl?dl=0)    |
| 10 million        |                                                                                                        |                                                                                                          | [0.5% 25kb](https://dl.dropboxusercontent.com/s/gsaxkgz0oh4ahgf/high-confidence.10million.25kb.pkl?dl=0) |                                                                                         |                                                                                               |
| 5 million         |                                                                                                        |                                                                                                          | [0.25% 25kb](https://dl.dropboxusercontent.com/s/10fbe85lpfabupw/high-confidence.5million.25kb.pkl?dl=0) |                                                                                         |                                                                                               |

To make it clear, let's download another Hi-C dataset from 4DN: https://data.4dnucleome.org/files-processed/4DNFI5IHU27G/@@download/4DNFI5IHU27G.mcool. Peakachu provides a handy function `peakachu depth` to extract the total number of intra-chromosomal pairs from *[cool URI](https://cooler.readthedocs.io/en/latest/concepts.html#uri-string)*:


```bash
peakachu depth -p 4DNFI5IHU27G.mcool:resolutions/1000000
```

 592991890

According to the table, for ~600 million intra-reads, we recommend using the 30% models in your prediction.

```bash
peakachu score_genome -r 10000 --balance -p 4DNFI5IHU27G.mcool:resolutions/10000 -O 4DNFI5IHU27G-peakachu-10kb-scores.bedpe -m high-confidence.600million.10kb.pkl
peakachu pool -r 10000 -i 4DNFI5IHU27G-peakachu-10kb-scores.bedpe -o 4DNFI5IHU27G-peakachu-10kb-loops.0.95.bedpe -t 0.95
```

# Not just Hi-C
Peakachu has been tested on Hi-C, Micrco-C, and DNA SPRITE contact maps with good results. For training sets, ChIA-PET, HiChIP, and PLAC-Seq have been tested. The purpose of this software is ultimately to facilitate the interpretation of results from multiple types of experiments, and the user is encouraged to apply Peakachu's training framework to newer approaches as they become available.
