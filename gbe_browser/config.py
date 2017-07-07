import os
import scidbpy


# == =
# == = Load = ==
# == =

SCIDB_INSTANCE_NUM = 2
GBE_DATA_PATH = '/home/scidb/GlobalBiobankEngine/gbe_data'


# -- -
# -- - Load: QC - --
# -- -
QC_PATH = os.path.join(GBE_DATA_PATH, 'qc')
QC_FILES = (
    {'file': os.path.join(QC_PATH, 'UKBioBiLallfreqSNPexclude.dat'),
     'header': 1},
    {'file': os.path.join(QC_PATH, 'ukb_ukbl_low_concordance.dat')}
)

QC_ARRAY = 'qc'


# -- -
# -- - Load: ICD - --
# -- -
ICD_GLOB = os.path.join(
    GBE_DATA_PATH, 'icdassoc', 'hybrid', 'c*.hybrid.rewrite.gz')
QT_GLOB = os.path.join(
    GBE_DATA_PATH, 'icdassoc', 'hybrid', 'c*.linear.rewrite.gz')

ICD_INFO_FILE = os.path.join(GBE_DATA_PATH, 'icdstats', 'icdinfo.txt')

ICD_INFO_ARRAY = 'icd_info'
ICD_INFO_SCHEMA = '<icd:string, Case:int64, Name:string>[icd_idx = 0:*:0:20]'

ICD_ARRAY = 'icd'
ICD_SCHEMA = """
  <icdind:      int64,
   affyid:      string,
   or_val:      double,
   se:          double,
   pvalue:      double,
   lor:         double,
   log10pvalue: double,
   l95or:       double,
   u95or:       double>
  [icd_idx   = 0:*:0:20;
   chrom     = 1:25:0:1;
   pos       = 0:*:0:10000000;
   pdecimal  = 0:3:0:1;
   synthetic = 0:999:0:1000]"""

ICD_PVALUE_MAP = dict(zip((.001, .0001, .00001), range(1, 4)))

# TODO 10K limit

ICD_INFO_STORE_QUERY = """
  store(
    redimension(
      apply(
        input({input_schema}, '{{fn}}', 0, 'CSV'),
        Case, int64(null),
        Name, string(null)),
      {icd_info_schema}),
    {icd_info_array})""".format(
        input_schema=ICD_INFO_SCHEMA.replace(', Case:int64, Name:string', ''),
        icd_info_schema=ICD_INFO_SCHEMA,
        icd_info_array=ICD_INFO_ARRAY)

ICD_INFO_INSERT_QUERY = """
  insert(
    join(
      project({icd_info_array}, icd),
      redimension(
        apply(
          index_lookup(
            aio_input('{path}', 'num_attributes=6'),
            project({icd_info_array}, icd),
            a0,
            icd_idx),
          Case, dcast(a1, int64(null)),
          Name, a2),
        {input_schema})),
    {icd_info_array})""".format(
        icd_info_array=ICD_INFO_ARRAY,
        input_schema=ICD_INFO_SCHEMA.replace('icd:string, ', ''),
        path=ICD_INFO_FILE)

ICD_INSERT_QUERY = """
  insert(
    redimension(
      apply(
        filter(
          index_lookup(
            aio_input(
              'paths={{paths}}',
              'instances={{instances}}',
              'num_attributes=12') as INPUT,
            {qc_array},
            INPUT.a2,
            is_in_filter),
          substr(a0, 0, 1) <> '#' and
          a6 = 'ADD' and
          is_in_filter is null and
          dcast(a9, double(null)) < .5 and
          a11 <> 'NA' and
          dcast(a11, double(null)) <> 0),
        icd_idx,     {{icd_idx_cond}},
        chrom,       int64(a0),
        pos,         int64(a1),
        pdecimal,    iif(dcast(a11, double(null)) < .00001, 3,
                      iif(dcast(a11, double(null)) < .0001, 2,
                       iif(dcast(a11, double(null)) < .001, 1, 0))),
        icdind,      int64(string(int64(a0) * 1e9 + int64(a1)) +
                           {{icdind_cond}}),
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
      {icd_array}),
    {icd_array})""".format(
        icd_array=ICD_ARRAY,
        qc_array=QC_ARRAY)

QT_INSERT_QUERY = """
  insert(
    redimension(
      apply(
        filter(
          index_lookup(
            aio_input(
              'paths={{paths}}',
              'instances={{instances}}',
              'num_attributes=12') as INPUT,
            {qc_array},
            INPUT.a2,
            is_in_filter),
          substr(a0, 0, 1) <> '#' and
          a5 = 'ADD' and
          is_in_filter is null and
          dcast(a8, double(null)) < .5 and
          a10 <> 'NA'),
        icd_idx,     {{icd_idx_cond}},
        chrom,       int64(a0),
        pos,         int64(a1),
        pdecimal,    iif(dcast(a10, double(null)) < .00001, 3,
                      iif(dcast(a10, double(null)) < .0001, 2,
                       iif(dcast(a10, double(null)) < .001, 1, 0))),
        icdind,      int64(string(int64(a0) * 1e9 + int64(a1)) +
                           {{icdind_cond}}),
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
      {icd_array}),
    {icd_array})""".format(
        icd_array=ICD_ARRAY,
        qc_array=QC_ARRAY)


# -- -
# -- - Load: GENE - --
# -- -
GENE_FILE = os.path.join(GBE_DATA_PATH, 'gencode.gtf.gz')

GENE_INDEX_ARRAY = 'gene_index'
GENE_INDEX_SCHEMA = '<gene_id:string>[gene_idx = 0:*:0:20]'

GENE_INDEX_STORE_QUERY = """
  store(
    redimension(
      apply(
        uniq(
          sort(
            project(
              apply(
                filter(
                  aio_input('{{path}}', 'num_attributes=9'),
                  substr(a0, 0, 1) <> '#'),
                gene_id, rsub(a8, 's/.*gene_id "([^.]*).*/$1/')),
              gene_id))),
        gene_idx, i),
      {gene_index_schema}),
    {gene_index_array})""".format(gene_index_array=GENE_INDEX_ARRAY,
                                  gene_index_schema=GENE_INDEX_SCHEMA)


GENE_ARRAY = 'gene'
GENE_SCHEMA = """
  <gene_name: string,
   strand:    string>
  [gene_idx  = 0:*:0:20;
   chrom     = 1:25:0:1;
   start     = 0:*:0:10000000;
   stop      = 0:*:0:10000000;
   synthetic = 0:199:0:200]"""

GENE_STORE_QUERY = """
  store(
    redimension(
      index_lookup(
        apply(
          filter(
            aio_input('{{path}}', 'num_attributes=9'),
            substr(a0, 0, 1) <> '#'),
          gene_id,      rsub(a8, 's/.*gene_id "([^.]*).*/$1/'),
          chrom,        iif(substr(a0, 3, 4) = 'X',
                            23,
                            iif(substr(a0, 3, 4) = 'Y',
                                24,
                                iif(substr(a0, 3, 4) = 'M',
                                    25,
                                    dcast(substr(a0, 3, 5), int64(null))))),
          start,        int64(a3) + 1,
          stop,         int64(a4) + 1,
          gene_name,    rsub(a8, 's/.*gene_name "([^"]*).*/$1/'),
          strand,       a6) as INPUT,
        {gene_index_array},
        INPUT.gene_id,
        gene_idx),
      {gene_schema}),
    {gene_array})""".format(gene_array=GENE_ARRAY,
                            gene_schema=GENE_SCHEMA,
                            gene_index_array=GENE_INDEX_ARRAY)


# -- -
# -- - Load: VARIANT - --
# -- -
VARIANT_FILE = os.path.join(GBE_DATA_PATH,
                            'icd10ukbb.ukbiobank.merge.sort.vcf.gz')

VARIANT_ARRAY = 'variant'
VARIANT_SCHEMA = """
  <rsid:         int64,
   ref:          string,
   alt:          string,
   site_quality: string,
   filter:       string,
   exac_nfe:     double,
   minpval:      double,
   minl10pval:   double,
   csq:          string>
  [chrom = 1:25:0:1;
   pos   = 0:*:0:10000000]"""

VARIANT_STORE_QUERY = """
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
        ref,          a3,
        alt,          a4,
        site_quality, a5,
        filter,       a6,
        exac_nfe,     dcast(rsub(a7, 's/.*EXAC_NFE=([^;]*).*/$1/'),
                            double(null)),
        minpval,      dcast(rsub(a7, 's/.*minpval=([^;]*).*/$1/'),
                            double(null)),
        minl10pval,   dcast(rsub(a7, 's/.*minl10pval=([^;]*).*/$1/'),
                            double(null)),
        csq,          rsub(a7, 's/.*CSQ=([^;]*).*/$1/')),
      {variant_array_schema}),
    {variant_array})""".format(variant_array=VARIANT_ARRAY,
                               variant_array_schema=VARIANT_SCHEMA)

VARIANT_GENE_ARRAY = 'variant_gene'
VARIANT_GENE_SCHEMA = """
  <val: int8 not null>
  [chrom    = 1:25:0:1;
   pos      = 0:*:0:10000000;
   gene_idx = 0:*:0:20]"""

VARIANT_GENE_STORE_QUERY = """
  store(
    redimension(
      index_lookup(
        apply(
          aio_input('{{path}}', 'num_attributes=3'),
          chrom,   int64(a0),
          pos,     int64(a1),
          gene_id, a2,
          val,     int8(0)) as INPUT,
        {gene_index_array},
        INPUT.gene_id,
        gene_idx),
      {variant_gene_schema}),
    {variant_gene_array})""".format(
        gene_index_array=GENE_INDEX_ARRAY,
        variant_gene_array=VARIANT_GENE_ARRAY,
        variant_gene_schema=VARIANT_GENE_SCHEMA)


# == =
# == = LOOKUP = ==
# == =

# -- -
# -- - Lookup: ICD - --
# -- -
ICD_LOOKUP_QUERY = """
  cross_join(
    {icd_array},
    filter({icd_info_array}, icd = '{{icd}}'),
    {icd_array}.icd_idx,
    {icd_info_array}.icd_idx)""".format(
        icd_array=ICD_ARRAY,
        icd_info_array=ICD_INFO_ARRAY)

ICD_PVALUE_LOOKUP_QUERY = """
  filter(
    cross_join(
      {icd_array},
      filter({icd_info_array}, icd = '{{icd}}'),
      {icd_array}.icd_idx,
      {icd_info_array}.icd_idx),
    pvalue < {{pvalue}})""".format(
        icd_array=ICD_ARRAY,
        icd_info_array=ICD_INFO_ARRAY)

ICD_CHROM_POS_LOOKUP_QUERY = """
  cross_join(
    between({icd_array}, null, {{chrom}}, {{pos}}, null, null,
                         null, {{chrom}}, {{pos}}, null, null),
    {icd_info_array},
    {icd_array}.icd_idx,
    {icd_info_array}.icd_idx)""".format(
        icd_array=ICD_ARRAY,
        icd_info_array=ICD_INFO_ARRAY)

ICD_X_INFO_SCHEMA = scidbpy.schema.Schema.fromstring(
    ICD_SCHEMA.replace(
        '>',
        ',{}>'.format(
            ICD_INFO_SCHEMA[ICD_INFO_SCHEMA.index('<') + 1:
                            ICD_INFO_SCHEMA.index('>')])))

ICD_VARIANT_LOOKUP_QUERY = """
  cross_join(
    project({variant_array},
            rsid,
            ref,
            alt,
            filter,
            exac_nfe,
            csq),
    cross_join(
        project(
          between({icd_array}, null, null, null, {{pdecimal}}, null,
                               null, null, null, null,         null),
          or_val,
          pvalue,
          log10pvalue),
        filter({icd_info_array}, icd = '{{icd}}'),
        {icd_array}.icd_idx,
        {icd_info_array}.icd_idx) as icd_join,
    {variant_array}.chrom,
    icd_join.chrom,
    {variant_array}.pos,
    icd_join.pos)""".format(
        icd_array=ICD_ARRAY,
        icd_info_array=ICD_INFO_ARRAY,
        variant_array=VARIANT_ARRAY)

VARIANT_X_ICD_X_INFO_SCHEMA = scidbpy.schema.Schema.fromstring("""
  <rsid:        int64,
   ref:         string,
   alt:         string,
   filter:      string,
   exac_nfe:    double,
   csq:         string,
   or_val:      double,
   pvalue:      double,
   log10pvalue: double,
   icd:         string,
   Case:        int64,
   Name:        string>
  [chrom     = 1:25:0:1;
   pos       = 0:*:0:10000000;
   icd_idx   = 0:*:0:20;
   pdecimal  = 0:3:0:1;
   synthetic = 0:999:0:1000]""")


# -- -
# -- - Lookup: GENE - --
# -- -

GENE_LOOKUP_QUERY = """
  cross_join(
    {gene_array},
    filter({gene_index_array}, gene_id = '{{gene_id}}'),
    {gene_array}.gene_idx,
    {gene_index_array}.gene_idx)""".format(
        gene_array=GENE_ARRAY,
        gene_index_array=GENE_INDEX_ARRAY)

GENE_LOOKUP_SCHEMA = scidbpy.schema.Schema.fromstring(
    GENE_SCHEMA.replace('>', ',gene_id:string>'))


# -- -
# -- - Lookup: VARIANT - --
# -- -
VARIANT_LOOKUP_QUERY = """
  between({variant_array}, {{chrom}}, {{pos}},
                     {{chrom}}, {{pos}})""".format(
    variant_array=VARIANT_ARRAY)

VARIANT_MULTI_LOOKUP_QUERY = """
  filter({variant_array}, {{chrom_pos_cond}})""".format(
    variant_array=VARIANT_ARRAY)

VARIANT_LOOKUP_SCHEMA = scidbpy.schema.Schema.fromstring(VARIANT_SCHEMA)

VARIANT_GENE_LOOKUP = """
  cross_join(
    {variant_array},
    cross_join(
      {variant_gene_array},
      filter({gene_index_array}, gene_id = '{{gene_id}}'),
      {variant_gene_array}.gene_idx,
      {gene_index_array}.gene_idx) as variant_gene_index,
    {variant_array}.chrom,
    variant_gene_index.chrom,
    {variant_array}.pos,
    variant_gene_index.pos)""".format(
        variant_array=VARIANT_ARRAY,
        variant_gene_array=VARIANT_GENE_ARRAY,
        gene_index_array=GENE_INDEX_ARRAY)

VARIANT_X_GENE_INDEX_SCHEMA = scidbpy.schema.Schema.fromstring(
    VARIANT_GENE_SCHEMA.replace(
        '<',
        '<{},'.format(
            VARIANT_SCHEMA[VARIANT_SCHEMA.index('<') + 1:
                           VARIANT_SCHEMA.index('>')]).replace(
                               '>',
                               ',{}>'.format(
                                   GENE_INDEX_SCHEMA[
                                       GENE_INDEX_SCHEMA.index('<') + 1:
                                       GENE_INDEX_SCHEMA.index('>')]))))

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
