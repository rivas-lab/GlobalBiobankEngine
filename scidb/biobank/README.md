# Loading Data

## Init SciDB

* SciDB Query: `init.afl`

<!-- -->

    iquery --afl --query-file init.afl

## Load Gene/Transcript/Exon

* Data:
  * `/opt/biobankengine/GlobalBioBankEngineRepo/gbe_data/canonical_transcripts.txt.gz`
  * `/opt/biobankengine/GlobalBioBankEngineRepo/gbe_data/omim_info.txt.gz`

<!-- -->

    gunzip --keep canonical_transcripts.txt.gz

    zcat omim_info.txt.gz                       \
    | tail --lines=+2                           \
    | cut --fields=1,3                          \
    | grep --invert --perl-regexp '\t$'         \
    | sort                                      \
    | uniq                                      \
    > omim_info-gene.txt

### HG19: gencode.gtf.gz

* Data: `/opt/biobankengine/GlobalBioBankEngineRepo/gbe_data/gencode.gtf.gz`

#### Gene: gencode-gene.gtf.gz

* SciDB Query: `load-gene.afl`

<!-- -->

    zcat gencode.gtf.gz                         \
    | grep --perl-regexp '\tgene\t'             \
    > gencode-gene.gtf

    iquery --afl --query-file load-gene.afl

#### Transcript: gencode-transcript.gtf

* SciDB Query: `load-transcript.afl`

<!-- -->

    zcat gencode.gtf.gz                         \
    | grep --perl-regexp '\ttranscript\t'       \
    > gencode-transcript.gtf

    iquery --afl --query-file load-transcript.afl

#### Exon: gencode-exon.gtf

* SciDB Query: `load-exon.afl`

<!-- -->

    zcat gencode.gtf.gz                                 \
    | tail --lines=+6                                   \
    | grep --invert --perl-regexp '\tgene\t'            \
    | grep --invert --perl-regexp '\ttranscript\t'      \
    > gencode-exon.gtf

    iquery --afl --query-file load-exon.afl

### HG38: gencode.v29.annotation.gtf.gz

* Data: `/opt/biobankengine/GlobalBioBankEngineRepo/gbe_data/gencode-gene.v29.annotation.gtf.gz`

#### Gene: gencode-gene.v29.annotation.gtf

* SciDB Query: `load-gene.afl`

<!-- -->

    zcat gencode.v29.annotation.gtf.gz          \
    | grep --perl-regexp '\tgene\t'             \
    > gencode-gene.v29.annotation.gtf

    iquery --afl --query-file load-gene.afl

#### Transcript: gencode-transcript.v29.annotation.gtf

* SciDB Query: `load-transcript.afl`

<!-- -->

    zcat gencode.gtf.gz                         \
    | grep --perl-regexp '\ttranscript\t'       \
    > gencode-transcript.gtf

    iquery --afl --query-file load-transcript.afl

#### Exon: gencode-exon.v29.annotation.gtf

* SciDB Query: `load-exon.afl`

<!-- -->

    cat gencode.v29.annotation.gtf.gz                   \
    | tail --lines=+6                                   \
    | grep --invert --perl-regexp '\tgene\t'            \
    | grep --invert --perl-regexp '\ttranscript\t'      \
    > gencode-exon.v29.annotation.gtf

    iquery --afl --query-file load-exon.afl

## Load Variants

### Variants Set: variant_filter_table.tsv.gz

* Data: `/opt/biobankengine/GlobalBioBankEngineRepo/gbe_data/qc/variant_filter_table.tsv.gz`
* SciDB Query: `load-variant-filter.afl`

<!-- -->

    zcat variant_filter_table.tsv.gz            \
    | tail --lines=+2                           \
    | sort --numeric-sort --key=1,1 --key=2,2   \
    | nl                                        \
        --line-increment=10                     \
        --number-format=ln                      \
        --number-width=1                        \
        --starting-line-number=0                \
    > variant_filter_table.scidb.tsv

    iquery --afl --query-file load-variant-filter.afl

### Variants Set: ukb_exm_spb-white_british-variant_annots_gbe.tsv.gz

* Data: `/opt/biobankengine/GlobalBioBankEngineRepo/gbe_data/ukb_exm_spb-white_british-variant_annots_gbe.tsv.gz`
* SciDB Query:
  * `create-variant-exome.afl`
  * `load-variant-exome.afl.tpl`

<!-- -->

    zcat ukb_exm_spb-white_british-variant_annots_gbe.tsv.gz     \
    | tail --lines=+2                                            \
    | sed 's/^X/23/'                                             \
    | sed 's/^Y/24/'                                             \
    | sort --numeric-sort --key=1,1 --key=2,2                    \
    | nl                                                         \
        --line-increment=10                                      \
        --number-format=ln                                       \
        --number-width=1                                         \
        --starting-line-number=0                                 \
    | split --numeric-suffixes=1 --lines=200000 -                \
        ukb_exm_spb-white_british-variant_annots_gbe.scidb.

    # 500MB files

    cat create-variant-exome.afl > load-variant-exome.afl
    for I in `seq --format='%02g' 53`; do        \
        I=$I envsubst                            \
        <  load-variant-exome.afl.tpl            \
        >> load-variant-exome.afl;               \
    done

    iquery --afl --query-file load-variant-exome.afl

    iquery --afl --query "remove_versions(LOAD.VARIANT, keep:1)"

### Variants Set: bbj_variant_annots.tsv.gz

* Data: `/opt/biobankengine/GlobalBioBankEngineRepo/gbe_data/bbj_variant_annots.tsv.gz`
* SciDB Query:
  * `create-variant-bbj.afl`
  * `load-variant-bbj.afl.tpl`

<!-- -->

    zcat bbj_variant_annots.tsv.gz    \
    | tail --lines=+2                                            \
    | sed 's/^X/23/'                                             \
    | sed 's/^Y/24/'                                             \
    | sort --numeric-sort --key=1,1 --key=2,2                    \
    | nl                                                         \
        --line-increment=10                                      \
        --number-format=ln                                       \
        --number-width=1                                         \
        --starting-line-number=0                                 \
    | split --numeric-suffixes=1 --lines=200000 -                \
        bbj_variant_annots_gbe.scidb.

    # 500MB files

    cat create-variant-bbj.afl > load-variant-bbj.afl
    for I in `seq --format='%02g' 77`; do        \
        I=$I envsubst                            \
        <  load-variant-bbj.afl.tpl            \
        >> load-variant-bbj.afl;               \
    done

    iquery --afl --query-file load-variant-bbj.afl

    iquery --afl --query "remove_versions(LOAD.VARIANT, keep:1)"

## Update Variants

### Update RSIDs

* Data:
  * HG 19: `/opt/biobankengine/GlobalBioBankEngineRepo/gbe_data/dbsnp151_hg19.txt.gz`
  * HG 38: `/opt/biobankengine/GlobalBioBankEngineRepo/gbe_data/dbsnp151_hg38.txt.gz`
* Python Script: `split-dbsnp.py`
* SciDB Query: `load-variant-rsid.afl`

<!-- -->

    zcat dbsnp151_hg19.txt.gz                   \
    | ./split-dbsnp.py                          \
    > dbsnp151_hg19.split.txt

    iquery --afl --query-file load-variant-rsid.afl

If the `VARIANT_UPDATED` array looks good, remove the `VARIANT` array
and rename `VARIANT_UPDATED` to `VARIANT`. Update the RSID index as
follows:

```R
$ R
> source('/usr/local/src/biobank/pkg/R/biobank.R')
> bb <- connect(username = "scidbadmin", password = "Paradigm4")
> bb$populate_rsid_index('RIVAS_HG19')
```

### Add Attributes

New attributes can be added to the `VARIANT` array as follows:

1. Prepare a `TAB` separated file with the new attributes. The first
   four columns of the file should be the chromosome, position,
   reference, and alternative for each updated variant.
1. Open and edit `update-variant.afl`
   1. Replace `NAMESPACE` with the namespace of the variant to update
   1. Replace `FILE_PATH` with the full path of the file with the new
      attributes
   1. Set `num_attributes` to the number of `TAB` separated columns in
      the file. At a minimum there should be four columns in the file
      to identify each variant
   1. List all the new attributes (besides the four required ones) and
      their position in the file (e.g., `a0`, `a1`, etc.). Perform any
      required type conversions
   1. List all the new attributes (besides the four required ones) and
      their type
1. Run the `update-variant.afl` file:

        iquery --afl --query-file update-variant.afl

   This results in a `VARIANT_UPDATED` array in the namespace of interest
1. Verify that the `VARIANT_UPDATED` has all the variants and the new
   fields:

        iquery --afl
        AFL% set_namespace(NAMESPACE);
        AFL% op_count(VARIANT);
        AFL% op_count(VARIANT_UPDATED);
        AFL% show(VARIANT_UPDATED);
        AFL% limit(VARIANT_UPDATED, 1);

1. It it checks out, remove the old `VARIANT` array and rename
   `VARIANT_UPDATED` to `VARIANT`:

        iquery --afl
        AFL% set_namespace(NAMESPACE);
        AFL% remove(VARIANT);
        AFL% rename(VARIANT_UPDATED, VARIANT);

1. Re-populate the variant fields using the following command. Replace
   `NAMESPACE` with the namespace of interest:

        NAMESPACE=NAMESPACE envsubst            \
        < populate-variant-fields.afl.tpl       \
        | iquery --afl

## Load Associations

### Association Set: data052019

* Data: `/opt/biobankengine/GlobalBioBankEngineRepo/gbe_data/data052019`
* SciDB R Script: `assoc.hg19.R`

<!-- -->

    Rscript assoc.hg19.R

Run with `nohup`:

    nohup Rscript assoc.hg19.R &

Build skip list from `nohup.out` log:

    grep "Done" nohup.out | cut --delimiter=" " --fields=3 > skip.hg19.txt

### Association Set: exome-20190514

* Data: `/opt/biobankengine/GlobalBioBankEngineRepo/gbe_data/exome-20190514`
* SciDB R Script: `assoc.hg38.R`

<!-- -->

    Rscript assoc.hg38.R

Run with `nohup`:

    nohup Rscript assoc.hg38.R &

Build skip list from `nohup.out` log:

    grep "Done" nohup.out | cut --delimiter=" " --fields=3 > skip.hg38.txt

### Association Set: bbj-20190822

* Data: `/opt/biobankengine/GlobalBioBankEngineRepo/gbe_data/bbj-20190822`
* SciDB R Script: `assoc.bbjhg19.R`

<!-- -->

    Rscript assoc.bbjhg19.R

Run with `nohup`:

    nohup Rscript assoc.bbjhg19.R &

Build skip list from `nohup.out` log:

    grep "Done" nohup.out | cut --delimiter=" " --fields=3 > skip.bbjhg19.txt 
    

## Update Phenotypes

* Data Directory `/opt/biobankengine/GlobalBioBankEngineRepo/gbe_data/`
  * Array: `icdinfo.shortname.tsv`
  * Exome: `icdinfo.shortnames.exome.tsv`
* SciDB Query: `load-phenotype.afl.tpl`

<!-- -->

    NAMESPACE=RIVAS_HG19                                \
    ASSOCIATION=ARRAY                                   \
    FILE=icdinfo.shortname.tsv                          \
    DATE=`date`                                         \
    envsubst < load-phenotype.afl.tpl | iquery --afl

If the `ASSOC_${ASSOCIATION}_FIELD_UPDATED` array looks good, remove
the `ASSOC_${ASSOCIATION}_FIELD` array and rename
`ASSOC_${ASSOCIATION}_FIELD_UPDATED` to `ASSOC_${ASSOCIATION}_FIELD`.

## Load Data in a New Namespace

1. Initialize the new namespace with the following command. Replace
   `NEW_NAMESPACE` with the name of the new namespace:

        NAMESPACE=NEW_NAMESPACE envsubst        \
        < init-new-namespace.afl.tpl            \
        | iquery --afl

1. Load gene/transcript/exon data using the earlier
   [instructions](#load-genetranscriptexon).

1. By default the gene/transcript/exon loaders will create new arrays
   in the `LOAD` namespace. These arrays have to be moved to the new
   namespace. Replace `NEW_NAMESPACE` with the name of the new
   namespace in the SciDB queries below:

        move_array_to_namespace(LOAD.GENE, NEW_NAMESPACE);
        move_array_to_namespace(LOAD.TRANSCRIPT, NEW_NAMESPACE);
        move_array_to_namespace(LOAD.EXON, NEW_NAMESPACE);

1. Load the variants in SciDB. See the earlier
   [examples](#load-variants).

1. By default the variant loaders will create a new `VARIANT` array in
   the `LOAD` namespace. This array has to be moved to the new
   namespace. Replace `NEW_NAMESPACE` with the name of the new
   namespace in the SciDB query below:

        move_array_to_namespace(LOAD.VARIANT, NEW_NAMESPACE);

1. Populate the variant fields using the following command. Replace
   `NEW_NAMESPACE` with the name of the new namespace:

        NAMESPACE=NEW_NAMESPACE envsubst        \
        < populate-variant-fields.afl.tpl       \
        | iquery --afl

1. Load the association data in SciDB. See the earlier
   [examples](#load-associations). Make sure to update the `R` script
   to use the new namespace. Also provide the name of the association
   set in the `R` script.

1. Update the phenotypes array in SciDB. See the earlier
   [example](#update-phenotypes)

1. Once all the steps are completed, check that the new namespace
   contains the required arrays. Replace `NEW_NAMESPACE` with the name
   of the new namespace in the SciDB queries below. In the
   output,`NEW_ASSOCIATION_SET` should be the name of the new
   association set:

        set_namespace(NEW_NAMESPACE);
        project(list(), name);
        {No} name
        {0} 'ASSOCIATION_SET'
        {1} 'ASSOC_NEW_ASSOCIATION_SET'
        {2} 'ASSOC_NEW_ASSOCIATION_SET_FIELD'
        {3} 'ASSOC_NEW_ASSOCIATION_SET_FIELD_MAP'
        {4} 'ASSOC_NEW_ASSOCIATION_SET_VARIANT'
        {5} 'EXON'
        {6} 'GENE'
        {7} 'NAMESPACE_DESCRIPTION'
        {8} 'PHENOTYPE_SET'
        {9} 'TRANSCRIPT'
        {10} 'VARIANT'
        {11} 'VARIANT_FIELDS'
        {12} 'VARIANT_RSID_INDEX'
