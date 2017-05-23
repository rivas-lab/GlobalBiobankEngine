# Global Biobank Engine
A framework for exploring and sharing genetic data from large genetic datasets.

## Setup

#### 1. Install requirements

```
pip install -r gbe_browser/requirements.txt
Rscript gbe_browser/install_R_requirements.R
```

#### 2. Extract data
```
tar -xvzf gbe_data.tar.gz
```

#### 3. Setup database

```
mkdir database
```

Run MongoDB.  For installation see [their website](https://docs.mongodb.com/manual/installation/)

```
mongod --dbpath database
```

Leave MongoDB running and open a new terminal.  In the new window run the following commands to load the database.

```
cd gbe_browser
python manage.py load_base_coverage
python manage.py load_variants_file
python manage.py load_gene_models
python manage.py load_dbsnp_file
python manage.py load_icd_info_stats
python manage.py load_icd_stats
```

#### 4. Start the GBE server

```
python gbe.py
```

#### 5. Use GBE

Open a web browser and navigate to `http://127.0.0.1:5000`




