# browser
Browser repository to make data and results publicly available


1) mkdir database

2) Run mongod --dbpath database.  (See [their website](https://docs.mongodb.com/manual/installation/) for installation instructions)

3) Open a new terminal window

4) Upload data to database, go to gbe_browser/

   * `python manage.py load_base_coverage`
   * `python manage.py load_variants_file`
   * `python manage.py load_gene_models`
   * `python manage.py load_dbsnp_file`
   * `python manage.py load_icd_info_stats`
   * `python manage.py load_icd_stats`

5) Run gbe.py server
   * `python gbe.py`

6) Open web browser, go to `http://127.0.0.1:5000`


