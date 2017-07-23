# biobankenginedev.stanford.edu

## Pull from GitHub.com

1. Git status:

```bash
/opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser$ git status
On branch scidb
Your branch is up-to-date with 'origin/scidb'.

Changes not staged for commit:
  (use "git add <file>..." to update what will be committed)
  (use "git checkout -- <file>..." to discard changes in working directory)

	modified:   client_secrets.json
	modified:   config.py

Untracked files:
...
```

2. Git stash, pull, and pop

```bash
/opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser$ git stash
Saved working directory and index state WIP on scidb: 15fc939 Remove second region example from timing function
HEAD is now at 15fc939 Remove second region example from timing function
```

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

```bash
/opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser$ git stash pop
On branch scidb
Your branch is up-to-date with 'origin/scidb'.

Changes not staged for commit:
  (use "git add <file>..." to update what will be committed)
  (use "git checkout -- <file>..." to discard changes in working directory)

	modified:   client_secrets.json
	modified:   config.py

Untracked files:
...
```

## Restart Flask

```bash
/opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser$ docker-compose restart flask
Restarting biobankengine_flask_1 ... done
```

## Load Data

```bash
/opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser$ python loader.py
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

### Load Individual Arrays

Running `loader.py` drops all arrays and load all the data files. To load individual arrays, one can individual functions of the `Loader` class like this:

```bash
/opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser$python -c 'import loader; loader.Loader().store_dbsnp()'
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

Because of the dependencies between index (e.g., `gene_index`) and
data (e.g., `gene`) it is required to reload the data arrays when the
index arrays are reloaded.

For example, to reload the `icdinfo.txt` file, run:

```bash
$python -c 'import loader; loader.Loader().insert_icd_info()'
```

For example, to reload variants file, run:

```bash
$python -c 'import loader; l = loader.Loader(); l.store_variant(); l.store_variant_gene(); l.store_variant_transcript()'
```
