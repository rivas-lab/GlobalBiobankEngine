import os
import scidbpy


SCIDB_INSTANCE_NUM = 2
GBE_DATA_PATH = '/home/scidb/GlobalBiobankEngine/gbe_data'


# -- -
# -- - QC - --
# -- -
QC_PATH = os.path.join(GBE_DATA_PATH, 'qc')
QC_FILES = (
    {'file': os.path.join(QC_PATH, 'UKBioBiLallfreqSNPexclude.dat'),
     'header': 1},
    {'file': os.path.join(QC_PATH, 'ukb_ukbl_low_concordance.dat')}
)

QC_ARRAY = 'qc'


# -- -
# -- - ICD - --
# -- -
ICD_GLOB = os.path.join(
    GBE_DATA_PATH, 'icdassoc', 'hybrid', 'c*.hybrid.rewrite.gz')
QT_GLOB = os.path.join(
    GBE_DATA_PATH, 'icdassoc', 'hybrid', 'c*.linear.rewrite.gz')

ICD_INFO_FILE = os.path.join(GBE_DATA_PATH, 'icdstats', 'icdinfo.txt')


ICD_INDEX_ARRAY = 'icd_index'
ICD_INDEX_SCHEMA = '<icd:string>[icd_id]'

ICD_ARRAY = 'icd'
ICD_SCHEMA = """
  <icdind:      int64,
   xpos:        int64,
   affyid:      string,
   or_val:      double,
   se:          double,
   pvalue:      double,
   lor:         double,
   log10pvalue: double,
   l95or:       double,
   u95or:       double>
  [icd_id    = 0:*:0:20;
   chrom     = 1:25:0:1;
   pos       = 0:*:0:10000000;
   synthetic = 0:999:0:1000]"""

ICD_INFO_ARRAY = 'icd_info'
ICD_INFO_SCHEMA = '<Case:int64, Name:string>[icd_id=0:*:0:1000000]'

ICD_LOAD_QUERY = """
  insert(
    redimension(
      apply(
        filter(
          index_lookup(
            aio_input(
              'paths={{paths}}',
              'instances={{instances}}',
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
        icd_id,      {{icd_id_cond}},
        chrom,       int64(a0),
        pos,         int64(a1),
        icdind,      int64(string(int64(a0) * 1e9 + int64(a1)) +
                           {{icdind_cond}}),
        xpos,        int64(a0) * int64(1e9) + int64(a1),
        affyid,      a2,
        or_val,      dcast(a8,  double(null)),
        se,          dcast(a9,  double(null)),
        pvalue,      dcast(a11, double(null)),
        lor,         log(dcast(a8, double(null))),
        log10pvalue, -log10(dcast(a11, double(null))),
        l95or,       exp(log(dcast(a8, double(null))) -
                         1.96 * dcast(a9, double(null))),
        u95or,       exp(log(dcast(a8, double(null))) +
                         1.96 * dcast(a9, double(null)))),
      {icd}),
    {icd})""".format(
        icd=ICD_ARRAY,
        qc=QC_ARRAY)

QT_LOAD_QUERY = """
  insert(
    redimension(
      apply(
        filter(
          index_lookup(
            aio_input(
              'paths={{paths}}',
              'instances={{instances}}',
              'num_attributes=12') as INPUT,
            {qc},
            INPUT.a2,
            is_in_filter),
          substr(a0, 0, 1) <> '#' and
          a5 = 'ADD' and
          is_in_filter is null and
          dcast(a8, double(null)) < .5 and
          a10 <> 'NA'),
        icd_id,      {{icd_id_cond}},
        chrom,       int64(a0),
        pos,         int64(a1),
        icdind,      int64(string(int64(a0) * 1e9 + int64(a1)) +
                           {{icdind_cond}}),
        xpos,        int64(a0) * int64(1e9) + int64(a1),
        affyid,      a1,
        or_val,      dcast(a7,  double(null)),
        se,          dcast(a8,  double(null)),
        pvalue,      dcast(a10, double(null)),
        lor,         dcast(a7, double(null)),
        log10pvalue, -log10(dcast(a10, double(null))),
        l95or,       exp(log(dcast(a7, double(null)))
                         - 1.96 * dcast(a8, double(null))),
        u95or,       exp(log(dcast(a7, double(null)))
                         + 1.96 * dcast(a8, double(null)))),
      {icd}),
    {icd})""".format(
        icd=ICD_ARRAY,
        qc=QC_ARRAY)

ICD_INFO_LOAD_QUERY = """
  store(
    redimension(
      apply(
        index_lookup(
          aio_input('{path}', 'num_attributes=6'),
          {icd_index},
          a0,
          icd_id),
        Case,     dcast(a1, int64(null)),
        Name, a2),
      {icd_info}, false),
    {icd_info})""".format(
        icd_info=ICD_INFO_ARRAY,
        icd_index=ICD_INDEX_ARRAY,
        path=ICD_INFO_FILE)


ICD_LOOKUP_QUERY = """
  filter(
    cross_join(
      {icd},
      filter({icd_index}, icd = '{{icd_id}}'),
      {icd}.icd_id,
      {icd_index}.icd_id),
    pvalue < {{cutoff}})""".format(
        icd=ICD_ARRAY,
        icd_index=ICD_INDEX_ARRAY)

ICD_LOOKUP_SCHEMA_INST = scidbpy.schema.Schema.fromstring(
    ICD_SCHEMA.replace('>', ',icd:string>'))

ICD_INFO_LOOKUP_QUERY = """
  cross_join(
    {icd_info},
    filter({icd_index}, icd = '{{icd_id}}'),
    {icd_info}.icd_id,
    {icd_index}.icd_id)""".format(
        icd_info=ICD_INFO_ARRAY,
        icd_index=ICD_INDEX_ARRAY)

ICD_INFO_LOOKUP_SCHEMA_INST = scidbpy.schema.Schema.fromstring(
    ICD_INFO_SCHEMA.replace('>', ',icd:string>'))


# -- -
# -- - VARIANT - --
# -- -
VARIANT_FILE = os.path.join(GBE_DATA_PATH,
                            'icd10ukbb.ukbiobank.merge.sort.vcf.gz')

VARIANT_ARRAY = 'variant'
VARIANT_SCHEMA = """
  <rsid:         int64,
   xpos:         int64,
   ref:          string,
   alt:          string,
   site_quality: string,
   filter:       string,
   minpval:      double,
   minl10pval:   double,
   csq:          string>
  [chrom     = 1:25:0:1;
   pos       = 0:*:0:10000000;
   synthetic = 0:999:0:1000]"""

VARIANT_LOAD_QUERY = """
  store(
    redimension(
      apply(
        filter(
          aio_input('{{path}}', 'num_attributes=8'),
          substr(a0, 0, 1) <> '#'),
        chrom,        int64(a0),
        pos,          int64(a1),
        rsid,         iif(strlen(a2) > 1,
                          int64(substr(a2, 2, 100)),
                          int64(null)),
        xpos,         int64(a0) * int64(1e9) + int64(a1),
        ref,          a3,
        alt,          a4,
        site_quality, a5,
        filter,       a6,
        minpval,      dcast(rsub(a7, 's/.*minpval=([^;]*).*/$1/'),
                            double(null)),
        minl10pval,   dcast(rsub(a7, 's/.*minl10pval=([^;]*).*/$1/'),
                            double(null)),
        csq,          rsub(a7, 's/.*CSQ=([^;]*).*/$1/')),
      {variant}),
    {variant})""".format(variant=VARIANT_ARRAY)

VARIANT_LOOKUP_QUERY = "filter({variant}, xpos = {{xpos}})".format(
    variant=VARIANT_ARRAY)
VARIANT_LOOKUP_SCHEMA_INST = scidbpy.schema.Schema.fromstring(VARIANT_SCHEMA)

VARIANT_CSQ = ('Allele',
               'Consequence',
               'IMPACT',
               'SYMBOL',
               'Gene',
               'Feature_type',
               'Feature',
               'BIOTYPE',
               'EXON',
               'INTRON',
               'HGVSc',
               'HGVSp',
               'cDNA_position',
               'CDS_position',
               'Protein_position',
               'Amino_acids',
               'Codons',
               'Existing_variation',
               'DISTANCE',
               'ALLELE_NUM',
               'STRAND',
               'VARIANT_CLASS',
               'MINIMISED',
               'SYMBOL_SOURCE',
               'HGNC_ID',
               'CANONICAL',
               'TSL',
               'APPRIS',
               'CCDS',
               'ENSP',
               'SWISSPROT',
               'TREMBL',
               'UNIPARC',
               'SIFT2',
               'SIFT',
               'PolyPhen',
               'DOMAINS',
               'HGVS_OFFSET',
               'GMAF',
               'AFR_MAF',
               'AMR_MAF',
               'EAS_MAF',
               'EUR_MAF',
               'SAS_MAF',
               'AA_MAF',
               'EA_MAF',
               'ExAC_MAF',
               'ExAC_Adj_MAF',
               'ExAC_AFR_MAF',
               'ExAC_AMR_MAF',
               'ExAC_EAS_MAF',
               'ExAC_FIN_MAF',
               'ExAC_NFE_MAF',
               'ExAC_OTH_MAF',
               'ExAC_SAS_MAF',
               'CLIN_SIG',
               'SOMATIC',
               'PHENO',
               'PUBMED',
               'MOTIF_NAME',
               'MOTIF_POS',
               'HIGH_INF_POS',
               'MOTIF_SCORE_CHANGE',
               'LoF',
               'LoF_filter',
               'LoF_flags',
               'LoF_info')
