# How to download files

## UK Biobank

At the moment of this writing, SNP microarray data is downloaded using the `gfetch` tool. The documentation here corresponds to [this page](https://biobank.ctsu.ox.ac.uk/crystal/ukb/docs/instruct_gfetch.html).

The tool `gfetch ` can be downloaded by running:

```
wget  -nd  biobank.ndph.ox.ac.uk/ukb/util/gfetch
```

### Genotyped calls

For the genotype data (`bed` file):
```
gfetch 22418 -c<CHROMOSOME>
```

For the sample file (`fam` file):
```
gfetch 22418 -c<CHROMOSOME> -m
```

For the SNP data (`bim` file, see [Resource 1963](https://biobank.ndph.ox.ac.uk/ukb/refer.cgi?id=1963)):
```
wget  -nd  biobank.ndph.ox.ac.uk/ukb/ukb/auxdata/ukb_snp_bim.tar
```

### Imputed genotypes

To download the BGEN files (chromosomes 1 through 22):

```
for CHR in `seq 1 22`; do  
  nohup ./gfetch 22828 -c${CHR} -ak11350.key > genotypes/imputed/full/download_bgen_chr${CHR}.log &  
done
```

To download the index files:
```
wget  -nd  biobank.ctsu.ox.ac.uk/ukb/ukb/auxdata/ukb_imp_bgi.tgz
mv ukb_imp_bgi.tgz genotypes/imputed/full
tar -xzvf  genotypes/imputed/full/ukb_imp_bgi.tgz
```

For the sample files, the following strategy is recommended in order to avoid data duplication:
```
nohup ./gfetch 22828 -m -c1 -ak11350.key > genotypes/imputed/full/download_bgen_chr${CHR}.log
mv ukb22828_c1_b0_v3_s487213.sample genotypes/imputed/full

for CHR in `seq 2 22`; do  
  ln -s ukb22828_c1_b0_v3_s487213.sample genotypes/imputed/full/ukb22828_c${CHR}_b0_v3_s487213.sample
done
```


### Whole exome sequencing data
To download WES data:

```
gfetch 23156 -c<CHROMOSOME> -b<BLOCK>
```

To determine how many blocks there are available for each chromosome, download the file `pvcf_blocks.txt`:
```
wget -nd  biobank.ndph.ox.ac.uk/ukb/ukb/auxdata/pvcf_blocks.txt
```
(See [Resource 837](https://biobank.ndph.ox.ac.uk/ukb/refer.cgi?id=837))

To download the `bim` files corresponding to the above files:

```
wget  -nd  biobank.ctsu.ox.ac.uk/ukb/ukb/auxdata/UKBexomeOQFEbim.zip
```

(See [Resource 200](https://biobank.ndph.ox.ac.uk/ukb/refer.cgi?id=200))


### Other genetic datasets

#### Relatedness
```
gfetch rel
```
