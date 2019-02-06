
# Introduction
Hi-C matrices are comprised of read counts of pair-wise loci in a chromosome conformation capture experiment. peakachu is a command line tool for training and applying random forest models which classify Hi-C pixels as anchor loci for DNA loops. Models are trained via a binary classification schema where loci from ChIA-PET experiments are used as a positive class and random loci as a negative class. peakachu is designed to mimic the input, requirements, and syntax of pyHICCUPS and it is encouraged to install both peak callers in the same virtual environment.

## Requirements
peakachu runs on python3.6 and uses the sklearn, cooler, and scipy packages. It's recommended to use conda to create a new environment like so:


```bash
conda create -n peak_caller python=3.6 scikit-learn=0.20.2 numpy scipy pandas h5py
conda install -n peak_caller -c bioconda cooler # if not available, add bioconda channel to your miniconda
source activate peak_caller
pip install -i https://test.pypi.org/simple/ peakachu  
```

Peakachu should now be available within this environment. 

## Example commands


```bash
# train forests with Hi-C dataset and ChIA-PET peaklist
peakachu train -p gm12878.cool::10000 -b CTCF_ChIAPET.bed -O out_folder
```


```bash
# apply models to Hi-C data per-chromosome
for model in out_folder/*pkl; do 
peakachu score_chromosome -p gm12878.cool::10000 -m $i -O out_folder; done
```


```bash
# apply one model on all chromosomes
peakachu score_genome -p gm12878.cool::10000 -m $i -O out_folder; done
```


```bash
# pool enriched regions and call loops to stdout
for bedfile in out_folder/*bed; do 
peakachu pool -i $bedfile > ${bedfile}.loops; done
```

## Interpretation of results
By default, the called loops are chosen by retaining the pixels with the highest probability scores, up to 2000 per chromosome. The output format is a .bedpe file with the 7th column holding the interaction probability determined by the classifier. Further filtering can be done via awk:


```bash
cat out_folder/*loops > gm12878.loops
wc -l gm12878.loops
```

       36512 gm12878.loops



```bash
awk '{if ($7 > .9) print $0}' gm12878.loops | sort -u > gm12878.filtered.loops
wc -l gm12878.filtered.loops
```

       11329 gm12878.filtered.loops

