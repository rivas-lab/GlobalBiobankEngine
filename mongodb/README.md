# Global Biobank Engine
A framework for exploring and sharing genetic data from large genetic datasets.
## Pulling data

```
git lfs pull --all

```
## Setup

#### 1. Install requirements

```
pip install -r gbe_browser/requirements.txt

```
We recommend using Python 2.7.

```
python2.7 -m ensurepip --upgrade

```
Then, running

```
pip2.7 install --user -r requirements.txt 
```
if you don't have sudo access.

#### 2. Extract data
```
tar -xvzf gbe_data.tar.gz
```

#### 3. Setup database

```
mkdir database
```

Run MongoDB.  For installation see [their website](https://docs.mongodb.com/manual/installation/).

For example, if installing on MAC

```
brew tap mongodb/brew
```

To install MongoDB

```
brew install mongodb-community@5.0
```

To get MongoDB running, run

```
mongod --dbpath database
```

Leave MongoDB running and open a new terminal.  In the new window run the following commands to load the database.

```
cd gbe_browser
python manage_mongodb.py load_base_coverage
python manage_mongodb.py load_variants_file
python manage_mongodb.py load_gene_models
python manage_mongodb.py load_dbsnp_file
python manage_mongodb.py load_icd_info_stats
python manage_mongodb.py load_icd_stats
```

#### 4. Start the GBE server

```
python gbe_mongodb.py
```

#### 5. Use GBE

Open a web browser and navigate to `http://127.0.0.1:5000`

## GBE data files

`SITES_VCFS = icd10ukbb.ukbiobank.merge.sort.vcf.gz` variant/sites file annotation,

`GENCODE_GTF = gencode.gtf.gz` Gencode reference transcript file,
 
`CANONICAL_TRANSCRIPT_FILE = canonical_transcripts.txt.gz` canonical transcript file,
 
`OMIM_FILE = omim_info.txt.gz` OMIM mapping file, 

`BASE_COVERAGE_FILES = gbe_data/coverage/Panel2016.*.coverage.txt.gz` variant mapping to maximum -log10(P),
    
`ICD_INFO_FILE = gbe_data/icdstats/icdinfo.txt` ICD info file described [here](https://github.com/rivas-lab/GlobalBiobankEngine/tree/master/gbe_data/icdstats),
   
`ICD_STATS_FILES = gbe_data/icdassoc/hybrid/*c*.hybrid.rewritewna.gz` logistic regression PLINK output file formats, 
 
`QT_STATS_FILES = gbe_data/icdassoc/hybrid/c*.linear.rewrite.gz` linear regression PLINK output file formats,
 
`DBSNP_FILE = gbe_data/dbsnp150.txt.gz` dbSNP150 variant files.


## Sample output

For example output loading MongoDB data see

```
python2.7 manage_mongodb.py load_base_coverage
starting gbe.py0
localhost
27017
MongoClient('localhost', 27017)
Dropped db.base_coverage
Loading subset 7 of 8 total: 2 contigs from Panel2016.all.coverage.txt.gz
Loading subset 0 of 8 total: 3 contigs from Panel2016.all.coverage.txt.gz
Loading subset 3 of 8 total: 3 contigs from Panel2016.all.coverage.txt.gz
Loading subset 5 of 8 total: 3 contigs from Panel2016.all.coverage.txt.gz
Loading subset 2 of 8 total: 3 contigs from Panel2016.all.coverage.txt.gz
Loading subset 4 of 8 total: 3 contigs from Panel2016.all.coverage.txt.gz
Loading subset 6 of 8 total: 2 contigs from Panel2016.all.coverage.txt.gz
Loading subset 1 of 8 total: 3 contigs from Panel2016.all.coverage.txt.gz
Finished loading subset 6 from  Panel2016.all.coverage.txt.gz (37254 records)
Finished loading subset 5 from  Panel2016.all.coverage.txt.gz (70726 records)
Finished loading subset 7 from  Panel2016.all.coverage.txt.gz (80922 records)
Finished loading subset 4 from  Panel2016.all.coverage.txt.gz (84872 records)
Loaded 100000 records from subset 3 of 8 from Panel2016.all.coverage.txt.gz (5 seconds)
Loaded 100000 records from subset 0 of 8 from Panel2016.all.coverage.txt.gz (5 seconds)
Loaded 100000 records from subset 2 of 8 from Panel2016.all.coverage.txt.gz (5 seconds)
Loaded 100000 records from subset 1 of 8 from Panel2016.all.coverage.txt.gz (5 seconds)
Finished loading subset 1 from  Panel2016.all.coverage.txt.gz (105899 records)
Finished loading subset 2 from  Panel2016.all.coverage.txt.gz (120250 records)
Finished loading subset 0 from  Panel2016.all.coverage.txt.gz (138737 records)
Finished loading subset 3 from  Panel2016.all.coverage.txt.gz (141095 records)
```

For second step

```
python2.7 manage_mongodb.py load_variants_file
starting gbe.py0
localhost
27017
MongoClient('localhost', 27017)
Dropped db.variants
Loading subset 5 of 8 total: 3 contigs from icd10ukbb.ukbiobank.merge.sort.vcf.gz
Loading subset 3 of 8 total: 3 contigs from icd10ukbb.ukbiobank.merge.sort.vcf.gz
Loading subset 1 of 8 total: 3 contigs from icd10ukbb.ukbiobank.merge.sort.vcf.gz
Loading subset 7 of 8 total: 2 contigs from icd10ukbb.ukbiobank.merge.sort.vcf.gz
Loading subset 4 of 8 total: 3 contigs from icd10ukbb.ukbiobank.merge.sort.vcf.gz
Loading subset 0 of 8 total: 3 contigs from icd10ukbb.ukbiobank.merge.sort.vcf.gz
Loading subset 2 of 8 total: 3 contigs from icd10ukbb.ukbiobank.merge.sort.vcf.gz
Loading subset 6 of 8 total: 2 contigs from icd10ukbb.ukbiobank.merge.sort.vcf.gz
Finished loading subset 6 from  icd10ukbb.ukbiobank.merge.sort.vcf.gz (37254 records)
Finished loading subset 5 from  icd10ukbb.ukbiobank.merge.sort.vcf.gz (70726 records)
Finished loading subset 4 from  icd10ukbb.ukbiobank.merge.sort.vcf.gz (84872 records)
Loaded 100000 records from subset 1 of 8 from icd10ukbb.ukbiobank.merge.sort.vcf.gz (42 seconds)
Finished loading subset 7 from  icd10ukbb.ukbiobank.merge.sort.vcf.gz (80922 records)
Finished loading subset 1 from  icd10ukbb.ukbiobank.merge.sort.vcf.gz (105899 records)
Loaded 100000 records from subset 3 of 8 from icd10ukbb.ukbiobank.merge.sort.vcf.gz (45 seconds)
Loaded 100000 records from subset 2 of 8 from icd10ukbb.ukbiobank.merge.sort.vcf.gz (50 seconds)
Loaded 100000 records from subset 0 of 8 from icd10ukbb.ukbiobank.merge.sort.vcf.gz (50 seconds)
Finished loading subset 2 from  icd10ukbb.ukbiobank.merge.sort.vcf.gz (120250 records)
Finished loading subset 3 from  icd10ukbb.ukbiobank.merge.sort.vcf.gz (141095 records)
Finished loading subset 0 from  icd10ukbb.ukbiobank.merge.sort.vcf.gz (138737 records)
```

For third step

```
python2.7 manage_mongodb.py load_gene_models
starting gbe.py0
localhost
27017
MongoClient('localhost', 27017)
Dropped db.genes, db.transcripts, and db.exons.
Done loading metadata. Took 1 seconds
Done loading genes. Took 13 seconds
Done indexing gene table. Took 2 seconds
Done loading transcripts. Took 12 seconds
Done indexing transcript table. Took 0 seconds
Done loading exons. Took 7 seconds
Done indexing exon table. Took 0 seconds

```

For fourth step

```
python2.7 manage_mongodb.py load_dbsnp_file
starting gbe.py0
localhost
27017
MongoClient('localhost', 27017)
Loading subset 5 of 8 total: 3 contigs from dbsnp150.txt.gz
Loading subset 7 of 8 total: 3 contigs from dbsnp150.txt.gz
Loading subset 6 of 8 total: 3 contigs from dbsnp150.txt.gz
Loading subset 0 of 8 total: 3 contigs from dbsnp150.txt.gz
Loading subset 4 of 8 total: 3 contigs from dbsnp150.txt.gz
Loading subset 3 of 8 total: 3 contigs from dbsnp150.txt.gz
Loading subset 1 of 8 total: 3 contigs from dbsnp150.txt.gz
Loading subset 2 of 8 total: 3 contigs from dbsnp150.txt.gz
```

For fifth step

Before running fifth step 
```
cp icdinfo.txt icdinfo.txt.bk 
cp icdinfo2018.txt icdinfo.txt
```

```
python2.7 manage_mongodb.py load_icd_info_stats
starting gbe.py0
localhost
27017
MongoClient('localhost', 27017)
Dropped db.icd_info
Mood_swings
Miserableness
Irritability
Sensitivity_/_hurt_feelings
Fed-up_feelings
Nervous_feelings
Worrier_/_anxious_feelings
Tense_/_'highly_strung'
Worry_too_long_after_embarrassment
Suffer_from_'nerves'
Loneliness,_isolation
```

For sixth step

```
...
f.PHENO1_c20.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c22.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c3.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c16.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c18.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c15.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c21.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c2.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c17.glm.logistic.hybrid.rewritewna.gz (100616 records)
Finished loading subset 4 from  ukb16698_v2.HC224.phe.all.ff.PHENO1_c9.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c4.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c11.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c13.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c6.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c8.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c5.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c10.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c12.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c7.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c19.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c14.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c1.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c20.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c22.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c3.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c16.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c18.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c15.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c21.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c2.glm.logistic.hybrid.rewritewna.gz, ukb16698_v2.HC224.phe.all.ff.PHENO1_c17.glm.logistic.hybrid.rewritewna.gz (155904 records)
```

