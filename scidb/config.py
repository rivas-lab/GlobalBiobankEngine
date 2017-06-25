import os.path


SCIDB_INSTANCE_NUM = 8
GBE_DATA_PATH = '/home/scidb/GlobalBiobankEngine/gbe_data'


# -- -
# -- - QC - --
# -- -
QC_PATH = os.path.join(GBE_DATA_PATH, 'qc')
QC_FILES = (
    {'filename': os.path.join(QC_PATH, 'UKBioBiLallfreqSNPexclude.dat'),
     'header': 1},
    {'filename': os.path.join(QC_PATH, 'ukb_ukbl_low_concordance.dat')}
)

QC_ARRAY = 'qc'


# -- -
# -- - ICD - --
# -- -
ICD_GLOB = os.path.join(GBE_DATA_PATH, 'icdassoc', 'hybrid', 'c*.hybrid.gz')

ICD_INDEX_ARRAY = 'icd_index'
ICD_INDEX_SCHEMA = "<icd:string>[icd_id]"

ICD_ARRAY = 'icd'
ICD_SCHEMA = """
  <affyid:      string,
   or_val:      double,
   se:          double,
   pvalue:      double,
   lor:         double,
   log10pvalue: double,
   l95or:       double,
   u95or:       double>
  [icd_id    = 0:*,20,0,
   chrom     = 1:25,1,0,
   pos       = 0:*,10000000,0,
   synthetic = 0:999,1000,0]"""
ICD_QUERY = """
  insert(
    redimension(
      apply(
        filter(
          index_lookup(
            aio_input(
              'paths={paths}',
              'instances={instances}',
              'num_attributes=12') as INPUT,
            {qc},
            INPUT.a2,
            is_in_filter),
          substr(a0, 0, 1) <> '#' and
          a6 = 'ADD' and
          is_in_filter is null and
          dcast(a9, double(null)) < .5 and
          a11 <> 'NA' and
          dcast(a11, double(null)) <> 0),
        icd_id,      {icd_id_cond},
        chrom,       int64(a0),
        pos,         int64(a1),
        affyid,      a2,
        or_val,      dcast(a8,  double(null)),
        se,          dcast(a9,  double(null)),
        pvalue,      dcast(a11, double(null)),
        lor,         log(dcast(a8, double(null))),
        log10pvalue, log10(dcast(a11, double(null))) * -1,
        l95or,       exp(log(dcast(a8, double(null))) -
                         1.96 * dcast(a9, double(null))),
        u95or,       exp(log(dcast(a8, double(null))) +
                         1.96 * dcast(a9, double(null)))),
      {icd}),
    {icd})"""
