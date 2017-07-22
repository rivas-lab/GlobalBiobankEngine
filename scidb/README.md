# biobankenginedev.stanford.edu

## Pull from GitHub.com

1. Git status:

```bash
rvernica@biobankenginedev:/opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser$ git status
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
rvernica@biobankenginedev:/opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser$ git stash
Saved working directory and index state WIP on scidb: 15fc939 Remove second region example from timing function
HEAD is now at 15fc939 Remove second region example from timing function
```

```bash
rvernica@biobankenginedev:/opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser$ git pull
remote: Counting objects: 5, done.
remote: Compressing objects: 100% (3/3), done.
remote: Total 5 (delta 1), reused 5 (delta 1), pack-reused 0
Unpacking objects: 100% (5/5), done.
From github.com:rivas-lab/GlobalBiobankEngine
   e95f1ad..73affc5  master     -> origin/master
Already up-to-date.
```

```bash
rvernica@biobankenginedev:/opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser$ git stash pop
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
rvernica@biobankenginedev:/opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser$ docker-compose restart flask
Restarting biobankengine_flask_1 ... done
```

## Load Data

```bash
rvernica@biobankenginedev:/opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser$ python loader.py
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
