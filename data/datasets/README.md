## UK Biobank

### How to download files

At the moment of this writing, SNP microarray data is downloaded using the `gfetch` tool. The documentation here corresponds to [this page](https://biobank.ctsu.ox.ac.uk/crystal/ukb/docs/instruct_gfetch.html).

The tool `gfetch ` can be downloaded by running:

```
wget  -nd  biobank.ndph.ox.ac.uk/ukb/util/gfetch
```

#### Genotyped calls

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

#### Imputed genotypes


#### Whole exome sequencing data
To download WES data:

```
gfetch 23156 -c<CHROMOSOME> -b<BLOCK>
```

To determine how many blocks there are available for each chromosome, download the file `pvcf_blocks.txt`:
```
wget  -nd  biobank.ndph.ox.ac.uk/ukb/ukb/auxdata/pvcf_blocks.txt
```
 
(See [Resource 837](https://biobank.ndph.ox.ac.uk/ukb/refer.cgi?id=837))

#### Relatedness


```
gfetch rel
```
