# biobankenginedev.stanford.edu

## Pull from GitHub.com

1. Check that branch is `scidb`:

```bash
/opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser$ git status
On branch scidb
Your branch is up-to-date with 'origin/scidb'.
...
```

2. Pull from GitHub.com

```bash
/opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser$ git pull
remote: Counting objects: 5, done.
remote: Compressing objects: 100% (3/3), done.
remote: Total 5 (delta 1), reused 5 (delta 1), pack-reused 0
Unpacking objects: 100% (5/5), done.
From github.com:rivas-lab/GlobalBiobankEngine
   e95f1ad..73affc5  master     -> origin/master
Already up-to-date.
```

## Restart Flask

```bash
/opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser$ docker-compose restart flask
Restarting biobankengine_flask_1 ... done
```

## Load Data

Usage:

```bash
/opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser$ python loader.py -h
usage: loader.py [-h] {all,icd-incremental}

SciDB Loader

positional arguments:
  {all,icd-incremental}
                        Loading choice

optional arguments:
  -h, --help            show this help message and exit
```

### All

Using the `all` argument, drops all arrays and reloads all the files:

```bash
/opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser$ python loader.py all
_new_conn:Starting new HTTP connection (1): localhost
_new_conn:Starting new HTTP connection (1): localhost
_new_conn:Starting new HTTP connection (1): localhost
_new_conn:Starting new HTTP connection (1): localhost
_new_conn:Starting new HTTP connection (1): localhost
_new_conn:Starting new HTTP connection (1): localhost
make_fifo:FIFO:/tmp/tmpPiBRXa/fifo
make_fifo:FIFO:/tmp/tmpVU4v3z/fifo
make_fifo:FIFO:/tmp/tmpwuBYHN/fifo
make_fifo:FIFO:/tmp/tmpEjEvnb/fifo
make_fifo:FIFO:/tmp/tmpMIBk_A/fifo
make_fifo:FIFO:/tmp/tmpqywPu0/fifo
Remove and recreate arrays, confirm with "y":
```

### ICD Incremental Load

Using the `incremental-icd` argument, the script looks for new ICDs
and loads them. If any new ICDs are found, they are listed in the
output. The `icdinfo.txt` file is reloaded if new ICSs are found:

```bash
/opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser$ python loader.py icd-incremental
...
set_icd_qt_lists:New ICDs found: cancer1075,cancer1077,HC265,HC189,HC184,HC186,HC180,HC182,HC286,HC328,HC121,HC282,HC280,HC281,HC320,FH1113,HC322,HC325,HC288,HC15,HC37,HC12,HC19,RH56,cancer1026,HC315,HC198,HC196,HC214,HC210,HC190
set_icd_qt_lists:New QTs found: INI94,INI95,INI1279
append_icd_info:Array:icd_info
insert_icd_info:Array:icd_info
...
```

If no new ICDs are found, a message is displayed and the script exits:

```
/opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser$ python loader.py icd-incremental
...
set_icd_qt_lists:New ICDs found:
set_icd_qt_lists:New QTs found:
Traceback (most recent call last):
  File "loader.py", line 568, in <module>
    loader.set_icd_qt_lists()
  File "loader.py", line 108, in set_icd_qt_lists
    config.QT_GLOB))
Exception: No new ICD or QT files found. Patterns used:
/opt/biobankengine/GlobalBioBankEngineRepo/gbe_data/icdassoc/hybrid/*c*.hybrid.rewritewna.gz
/opt/biobankengine/GlobalBioBankEngineRepo/gbe_data/icdassoc/hybrid/*c*.linear.rewritewna.gz
```

### Load Individual Arrays

To load individual arrays, one can call individual functions of the
`Loader` class like this:

```bash
/opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser$ python -c 'import loader; loader.Loader().store_dbsnp()'
make_fifo:FIFO:/tmp/tmp4FEOze/fifo
make_fifo:FIFO:/tmp/tmp7TaW3p/fifo
make_pipe:Spawn:zcat /opt/biobankengine/GlobalBioBankEngineRepo/gbe_data/dbsnp150.txt.gz > /tmp/tmp4FEOze/fifo pid:3173
store_dbsnp:Query:running...
store_dbsnp:Query:done
store_dbsnp:Pipe:return code:0
store_dbsnp:Array:dbsnp_by_rsid
store_dbsnp:Query:running...
store_dbsnp:Query:done
store_dbsnp:Array:dbsnp_by_chrom_pos
remove_fifo:Remove:/tmp/tmp4FEOze/fifo
remove_fifo:Remove:/tmp/tmp7TaW3p/fifo
```

For the list of functions available see the `loader.py` file:

```bash
/opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser$ tail -25 loader.py
if __name__ == '__main__':
    loader = Loader()
    loader.remove_arrays()

    loader.store_qc()
    loader.set_icd_qt_lists()
    loader.append_icd_info()
    loader.insert_icd_info()
    loader.insert_icd()
    loader.insert_qt()

    loader.store_gene_index()
    loader.store_transcript_index()
    loader.store_dbnsfp()
    loader.store_canonical()
    loader.store_omim()
    loader.store_gene()
    loader.store_transcript()
    loader.store_exon()

    loader.store_variant()
    loader.store_variant_gene()
    loader.store_variant_transcript()

    loader.store_coverage()
    loader.store_dbsnp()
```

Because of the dependencies between index (e.g., `gene_index`) and
data (e.g., `gene`) arrays, it is required to reload the data arrays
when the index arrays are reloaded.

For example, to reload the `icdinfo.txt` file, run:

```bash
$ python -c 'import loader; loader.Loader().insert_icd_info()'
```

For example, to reload the variants file, run:

```bash
$ python -c 'import loader; l = loader.Loader(); l.store_variant(); l.store_variant_gene(); l.store_variant_transcript()'
```

### Update Input File Names

The input file names and file name patterns are located in
`config.py`. The variable names used end with `_FILE`, `_FILES`, or
`_GLOB`. For example the ICD file patterns are set in `ICD_GLOB` and
`QT_GLOB` variables, e.g.:

```python
ICD_GLOB = os.path.join(
    GBE_DATA_PATH, 'icdassoc', 'hybrid', '*c*.hybrid.rewritewna.gz')
QT_GLOB = os.path.join(
    GBE_DATA_PATH, 'icdassoc', 'hybrid', '*c*.linear.rewritewna.gz')
```
