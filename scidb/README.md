## Directory structure for GBE codebase

```
scidb/biobank/
```

directory containing instructions for initializing and loading data to SciDB.

```
scidb/streaming/
```

Sample codebase for streaming functions to the database. 

```
scidb/gbe_browser
```

Code for GBE web browser and html templates.

```
scidb/rivas/
```
Initial browser data loading instructions.

## GBE data files

`gbe_data/icdassoc/hybrid/` PLINK METAL results input files. 

`gbe_data/icdstats/icdinfo.txt` ICD info file. Format described [here](https://github.com/rivas-lab/GlobalBiobankEngine/tree/master/gbe_data/icdstats).

`gbe_data/gencode.gtf.gz` Gencode reference transcript and gene file.

`gbe_data/canonical_transcripts.txt.gz` Canonical transcripts files. 

`gbe_data/dbNSFP2.6_gene.gz` functional prediction of all nonsynonymous variants. 

`gbe_data/omim_info.txt.gz` OMIM reference file.

`gbe_data/icd10ukbb.ukbiobank.combined.vcf.gz` Variant annotation reference file. 

`gbe_data/coverage/Panel2016.all.coverage.txt.gz` maximum -log10(P) file across traits. Schema described [here](https://github.com/rivas-lab/GlobalBiobankEngine/blob/04771d9e9c9b978204607735d6a8dface6c7049a/scidb/gbe_browser/config.py#L751).

`gbe_data/dbsnp150.txt.gz` dbSNP 150 file.

