
# HiForest Introduction
Hi-C matrices are comprised of read counts of pair-wise loci in a chromosome conformation capture experiment. HiForest is a python command line tool for training and applying random forest models which classify Hi-C pixels as anchor loci for DNA loops. Models are trained via a binary classification schema where loci from ChIA-PET experiments are used as a positive class and random loci as a negative class. HiForest is designed to mimic the input, requirements, and syntax of pyHICCUPS and it is encouraged to install both peak callers in the same virtual environment.

## Requirements
HiForest runs on python3.6 and uses the sklearn, cooler, and scipy packages. It's recommended to use conda to create a new environment like so:


```bash
conda create -n peak_caller python=3.6 scikit-learn=0.20.2 numpy scipy pandas h5py
conda install -n peak_caller -c bioconda cooler
source activate peak_caller
pip install -i https://test.pypi.org/simple/ hiforest
```

## Example commands


```bash
# train forests with Hi-C dataset and ChIA-PET peaklist
hiforest train -p gm12878.cool::10000 -b CTCF_ChIAPET.bed -O out_folder
```


```bash
# apply models to Hi-C data per-chromosome
for model in out_folder/*pkl; do 
hiforest score_chromosome -p gm12878.cool::10000 -m $i -O out_folder; done
```


```bash
# apply one model on all chromosomes
hiforest score_genome -p gm12878.cool::10000 -m $i -O out_folder; done
```


```bash
# pool enriched regions and call loops
for bedfile in out_folder/*bed; do 
hiforest pool -i $bedfile > ${bedfile}.loops; done
```
