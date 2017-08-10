import os
import scidbpy


XOFF = int(1e9)

# == =
# == = Load = ==
# == =

SCIDB_INSTANCE_NUM = 6
GBE_DATA_PATH = '/opt/biobankengine/GlobalBioBankEngineRepo/gbe_data'


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
    GBE_DATA_PATH, 'icdassoc', 'hybrid', '*c*.hybrid.rewritewna.gz')
QT_GLOB = os.path.join(
    GBE_DATA_PATH, 'icdassoc', 'hybrid', '*c*.linear.rewritewna.gz')

ICD_INFO_FILE = os.path.join(GBE_DATA_PATH, 'icdstats', 'icdinfo.txt')

ICD_INFO_ARRAY = 'icd_info'
ICD_INFO_SCHEMA = '<icd:string, Case:int64, Name:string>[icd_idx = 0:*:0:80]'

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
  [icd_idx   = 0:*:0:80;
   chrom     = 1:25:0:1;
   pos       = 0:*:0:10000000;
   pdecimal  = 0:3:0:1;
   synthetic = 0:999:0:1000]"""

AFFYID_INDEX_ARRAY = 'affyid_index'
AFFYID_INDEX_SCHEMA = '<affyid:string>[affyid_idx = 0:*:0:100000]'

ICD_AFFYID_ARRAY = 'icd_affyid'
ICD_AFFYID_SCHEMA = """
  <chrom: int64,
   pos:   int64>
  [affyid_idx = 0:*:0:100000]"""

ICD_PVALUE_MAP = dict(zip((.001, .0001, .00001), range(1, 4)))

# TODO 10K limit

ICD_INFO_APPEND_QUERY = """
  insert(
    redimension(
      apply(
        input({input_schema}, '{{fn}}', 0, 'CSV'),
        Case, int64(null),
        Name, string(null)),
      {icd_info_schema}),
    {icd_info_array})""".format(
        input_schema=ICD_INFO_SCHEMA.replace(
            ', Case:int64, Name:string', '').replace(
                '0:*', '{start}:{stop}'),
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
        affyid,      a2,
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

AFFYID_INDEX_STORE_QUERY = """
  store(
    redimension(
      unpack(
        grouped_aggregate({icd_array}, count(*), affyid),
        affyid_idx),
      {affyid_index_schema}),
    {affyid_index_array})""".format(icd_array=ICD_ARRAY,
                                    affyid_index_array=AFFYID_INDEX_ARRAY,
                                    affyid_index_schema=AFFYID_INDEX_SCHEMA)

ICD_AFFYID_STORE_QUERY = """
  store(
    redimension(
      equi_join(
        project({icd_array}, affyid),
        {affyid_index_array},
        'left_names=affyid',
        'right_names=affyid',
        'keep_dimensions=1',
        'algorithm=hash_replicate_right'),
      {icd_affyid_schema},
      false),
    {icd_affyid_array})""".format(icd_array=ICD_ARRAY,
                                  affyid_index_array=AFFYID_INDEX_ARRAY,
                                  icd_affyid_schema=ICD_AFFYID_SCHEMA,
                                  icd_affyid_array=ICD_AFFYID_ARRAY)


# -- -
# -- - Load: GENE - --
# -- -
GENE_FILE = os.path.join(GBE_DATA_PATH, 'gencode.gtf.gz')
CANONICAL_FILE = os.path.join(GBE_DATA_PATH, 'canonical_transcripts.txt.gz')
DBNSFP_FILE = os.path.join(GBE_DATA_PATH, 'dbNSFP2.6_gene.gz')
OMIM_FILE = os.path.join(GBE_DATA_PATH, 'omim_info.txt.gz')

GENE_INDEX_ARRAY = 'gene_index'
GENE_INDEX_SCHEMA = '<gene_id: string>[gene_idx = 0:*:0:10000]'
GENE_INDEX_SCHEMA_OBJ = scidbpy.schema.Schema.fromstring(GENE_INDEX_SCHEMA)

GENE_INDEX_STORE_QUERY = """
  store(
    redimension(
      apply(
        aio_input('{{path}}', 'num_attributes=1'),
        gene_idx, tuple_no,
        gene_id,  rsub(a0, 's/"([^.]*).*/$1/')),
      {gene_index_schema}),
    {gene_index_array})""".format(gene_index_schema=GENE_INDEX_SCHEMA,
                                  gene_index_array=GENE_INDEX_ARRAY)

TRANSCRIPT_INDEX_ARRAY = 'transcript_index'
TRANSCRIPT_INDEX_SCHEMA = """
  <transcript_id:string>[transcript_idx = 0:*:0:10000]"""

TRANSCRIPT_INDEX_STORE_QUERY = """
  store(
    redimension(
      apply(
        aio_input('{{path}}', 'num_attributes=1'),
        transcript_idx, tuple_no,
        transcript_id,  rsub(a0, 's/"([^.]*).*/$1/')),
      {transcript_index_schema}),
    {transcript_index_array})""".format(
        transcript_index_schema=TRANSCRIPT_INDEX_SCHEMA,
        transcript_index_array=TRANSCRIPT_INDEX_ARRAY)

DBNSFP_ARRAY = 'dbnsfp'
DBNSFP_SCHEMA = '<full_gene_name: string>[gene_idx = 0:*:0:10000]'

DBNSFP_STORE_QUERY = """
  store(
    redimension(
      index_lookup(
        apply(
          filter(
            aio_input('{{path}}', 'num_attributes=14', 'header=1'),
            a1 <> '.' and a0 <> a12),
          full_gene_name, a12),
        {gene_index_array},
        a1,
        gene_idx),
      {dbnsfp_schema},
      false),
    {dbnsfp_array})""".format(gene_index_array=GENE_INDEX_ARRAY,
                              dbnsfp_schema=DBNSFP_SCHEMA,
                              dbnsfp_array=DBNSFP_ARRAY)

CANONICAL_ARRAY = 'canonical'
CANONICAL_SCHEMA = '<transcript_idx: int64>[gene_idx = 0:*:0:10000]'

CANONICAL_STORE_QUERY = """
  store(
    redimension(
      index_lookup(
        index_lookup(
          aio_input('{{path}}', 'num_attributes=2'),
          {gene_index_array},
          a0,
          gene_idx),
        {transcript_index_array},
        a1,
        transcript_idx),
      {canonical_schema},
      false),
    {canonical_array})""".format(gene_index_array=GENE_INDEX_ARRAY,
                                 transcript_index_array=TRANSCRIPT_INDEX_ARRAY,
                                 canonical_schema=CANONICAL_SCHEMA,
                                 canonical_array=CANONICAL_ARRAY)

OMIM_ARRAY = 'omim'
OMIM_SCHEMA = '<omim_accession: string>[gene_idx = 0:*:0:10000]'

OMIM_STORE_QUERY = """
  store(
    redimension(
      index_lookup(
        apply(
          filter(
            aio_input('{{path}}', 'num_attributes=4'),
            a2 <> ''),
          omim_accession, a2),
        {gene_index_array},
        a0,
        gene_idx),
      {omim_schema},
      false),
    {omim_array})""".format(gene_index_array=GENE_INDEX_ARRAY,
                            omim_schema=OMIM_SCHEMA,
                            omim_array=OMIM_ARRAY)

GENE_ARRAY = 'gene'
GENE_SCHEMA = """
  <gene_name:         string,
   strand:            string,
   full_gene_name:    string,
   omim_accession:    string,
   transcript_idx:    int64,
   chrom:             int64,
   start:             int64,
   stop:              int64,
   c_transcript_info: string,
   transcript_info:   string,
   exon_info:         string>
  [gene_idx       = 0:*:0:10000]"""

GENE_STORE_QUERY = """
  store(
    redimension(
      join(
        join(
          join(
            redimension(
              index_lookup(
                apply(
                  aio_input('{{path}}', 'num_attributes=9'),
                  g_id,      rsub(a8, 's/.*gene_id "([^.]*).*/$1/'),
                  chrom,     iif(substr(a0, 3, 4) = 'X',
                               23,
                               iif(substr(a0, 3, 4) = 'Y',
                                 24,
                                 iif(substr(a0, 3, 4) = 'M',
                                   25,
                                   dcast(substr(a0, 3, 5), int64(null))))),
                  start,     int64(a3),
                  stop,      int64(a4),
                  gene_name, rsub(a8, 's/.*gene_name "([^"]*).*/$1/'),
                  strand,    a6,
                  c_transcript_info, string(null),
                  transcript_info,   string(null),
                  exon_info,         string(null)),
                {gene_index_array},
                g_id,
                gene_idx),

              <gene_name: string,
               strand:    string,
               chrom:     int64,
               start:     int64,
               stop:      int64,
               c_transcript_info: string,
               transcript_info:   string,
               exon_info:         string>
              [gene_idx = 0:*:0:10000]),

            merge({dbnsfp_array},
                  project(apply({gene_index_array}, empty, string(null)),
                          empty))),
          merge({canonical_array},
                project(apply({gene_index_array}, empty, int64(null)),
                        empty))),
        merge({omim_array},
              project(apply({gene_index_array}, empty, string(null)),
                      empty))),
      {gene_schema}),
    {gene_array})""".format(gene_array=GENE_ARRAY,
                            gene_schema=GENE_SCHEMA,
                            gene_index_array=GENE_INDEX_ARRAY,
                            dbnsfp_array=DBNSFP_ARRAY,
                            canonical_array=CANONICAL_ARRAY,
                            omim_array=OMIM_ARRAY)

TRANSCRIPT_ARRAY = 'transcript'
TRANSCRIPT_SCHEMA = """
  <strand: string,
   chrom:  int64,
   start:  int64,
   stop:   int64>
  [gene_idx       = 0:*:0:10000;
   transcript_idx = 0:*:0:10000]"""

TRANSCRIPT_STORE_QUERY = """
  store(
    redimension(
      index_lookup(
        index_lookup(
          apply(
            aio_input('{{path}}', 'num_attributes=9'),
            g_id,   rsub(a8, 's/.*gene_id "([^.]*).*/$1/'),
            t_id,   rsub(a8, 's/.*transcript_id "([^.]*).*/$1/'),
            chrom,  iif(substr(a0, 3, 4) = 'X',
                        23,
                        iif(substr(a0, 3, 4) = 'Y',
                            24,
                            iif(substr(a0, 3, 4) = 'M',
                                25,
                                dcast(substr(a0, 3, 5), int64(null))))),
            start,  int64(a3),
            stop,   int64(a4),
            strand, a6),
          {gene_index_array},
          g_id,
          gene_idx),
        {transcript_index_array},
        t_id,
        transcript_idx),
      {transcript_schema}),
    {transcript_array})""".format(
        transcript_array=TRANSCRIPT_ARRAY,
        transcript_schema=TRANSCRIPT_SCHEMA,
        gene_index_array=GENE_INDEX_ARRAY,
        transcript_index_array=TRANSCRIPT_INDEX_ARRAY)

EXON_ARRAY = 'exon'
EXON_SCHEMA = """
  <feature_type: string>
  [gene_idx       = 0:*:0:20;
   transcript_idx = 0:*:0:20;
   chrom          = 1:25:0:1;
   start          = 0:*:0:10000000;
   stop           = 0:*:0:10000000;
   synthetic      = 0:199:0:200]"""
EXON_SCHEMA_OBJ = scidbpy.schema.Schema.fromstring(EXON_SCHEMA)

EXON_STORE_QUERY = """
  store(
    redimension(
      index_lookup(
        index_lookup(
          apply(
            aio_input('{{path}}', 'num_attributes=9'),
            g_id,         rsub(a8, 's/.*gene_id "([^.]*).*/$1/'),
            t_id,         rsub(a8, 's/.*transcript_id "([^.]*).*/$1/'),
            chrom,        iif(substr(a0, 3, 4) = 'X',
                              23,
                              iif(substr(a0, 3, 4) = 'Y',
                                  24,
                                  iif(substr(a0, 3, 4) = 'M',
                                      25,
                                      dcast(substr(a0, 3, 5), int64(null))))),
            start,        int64(a3),
            stop,         int64(a4),
            feature_type, a2),
          {gene_index_array},
          g_id,
          gene_idx),
        {transcript_index_array},
        t_id,
        transcript_idx),
      {exon_schema}),
    {exon_array})""".format(
        exon_array=EXON_ARRAY,
        exon_schema=EXON_SCHEMA,
        gene_index_array=GENE_INDEX_ARRAY,
        transcript_index_array=TRANSCRIPT_INDEX_ARRAY)

GENE_C_TRANSCRIPT_INFO_STORE_QUERY = """
  store(
    redimension(
      equi_join(
        project({gene_array},
                gene_name,
                strand,
                full_gene_name,
                omim_accession,
                transcript_idx,
                chrom,
                start,
                stop,
                transcript_info,
                exon_info),
        project(
          apply(
            equi_join(
              {transcript_index_array},
              equi_join(
                {transcript_array},
                apply(
                  project({gene_array}, transcript_idx),
                  g_idx,
                  gene_idx),
                'left_names=transcript_idx',
                'right_names=transcript_idx',
                'keep_dimensions=1',
                'algorithm=hash_replicate_right'),
              'left_names=transcript_idx',
              'right_names=transcript_idx',
              'algorithm=hash_replicate_right'),
            c_transcript_info,
            transcript_id + ':' + strand + ':' + string(chrom) + ':' +
            string(start) + ':' + string(stop)),
          g_idx,
          c_transcript_info),
        'left_names=gene_idx',
        'right_names=g_idx',
        'algorithm=hash_replicate_right'),
      {gene_array}),
    {gene_array})""".format(gene_array=GENE_ARRAY,
                            transcript_array=TRANSCRIPT_ARRAY,
                            transcript_index_array=TRANSCRIPT_INDEX_ARRAY)

GENE_TRANSCRIPT_INFO_STORE_QUERY = """
  store(
    redimension(
      equi_join(
        project({gene_array},
                gene_name,
                strand,
                full_gene_name,
                omim_accession,
                transcript_idx,
                chrom,
                start,
                stop,
                c_transcript_info,
                exon_info),
        aggregate(
          redimension(
            project(
              apply(
                equi_join(
                  {transcript_index_array},
                  equi_join(
                    {transcript_array},
                    project({gene_array}, gene_name),
                    'left_names=gene_idx',
                    'right_names=gene_idx',
                    'keep_dimensions=1',
                    'algorithm=hash_replicate_right'),
                  'left_names=transcript_idx',
                  'right_names=transcript_idx',
                  'algorithm=hash_replicate_right'),
                transcript_info,
                transcript_id + ':' + strand + ':' + string(chrom) + ':' +
                string(start) + ':' + string(stop) + ';'),
              gene_idx,
              transcript_info),
            <transcript_info: string>
            [gene_idx  = 0:*:0:1000000;
             synthetic = 0:*:0:1000]),
          sum(transcript_info) as transcript_info,
          gene_idx),
        'left_names=gene_idx',
        'right_names=gene_idx'),
      {gene_array}),
    {gene_array})""".format(gene_array=GENE_ARRAY,
                            transcript_array=TRANSCRIPT_ARRAY,
                            transcript_index_array=TRANSCRIPT_INDEX_ARRAY)

GENE_EXON_INFO_STORE_QUERY = """
  store(
    redimension(
      equi_join(
        project({gene_array},
                gene_name,
                strand,
                full_gene_name,
                omim_accession,
                transcript_idx,
                chrom,
                start,
                stop,
                c_transcript_info,
                transcript_info),
        aggregate(
          redimension(
            project(
              apply(
                equi_join(
                  {exon_array},
                  apply(
                    project({gene_array}, transcript_idx),
                    g_idx,
                    gene_idx),
                  'left_names=transcript_idx',
                  'right_names=transcript_idx',
                  'keep_dimensions=1',
                  'algorithm=hash_replicate_right'),
                exon_info,
                feature_type + ':' + string(chrom) + ':' +
                string(start) + ':' + string(stop) + ';'),
              g_idx,
              exon_info),
            <exon_info: string>
            [g_idx  = 0:*:0:1000000;
             synthetic = 0:*:0:1000]),
          sum(exon_info) as exon_info,
          g_idx),
        'left_names=gene_idx',
        'right_names=g_idx'),
      {gene_array}),
    {gene_array})""".format(gene_array=GENE_ARRAY,
                            exon_array=EXON_ARRAY)

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
   minicd:       string,
   minpval:      double,
   minor:        double,
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
        minicd,       rsub(a7, 's/.*minicd=([^;]*).*/$1/'),
        minpval,      dcast(rsub(a7, 's/.*minpval=([^;]*).*/$1/'),
                            double(null)),
        minor,        dcast(rsub(a7, 's/.*minor=([^;]*).*/$1/'),
                            double(null)),
        minl10pval,   dcast(rsub(a7, 's/.*minl10pval=([^;]*).*/$1/'),
                            double(null)),
        csq,          rsub(a7, 's/.*CSQ=([^;]*).*/$1/')),
      {variant_array_schema}),
    {variant_array})""".format(variant_array=VARIANT_ARRAY,
                               variant_array_schema=VARIANT_SCHEMA)

VARIANT_GENE_ARRAY = 'variant_gene'
VARIANT_GENE_SCHEMA = """
  <notused: int8 not null>
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
          notused, int8(0)),
        {gene_index_array},
        a2,
        gene_idx),
      {variant_gene_schema}),
    {variant_gene_array})""".format(
        gene_index_array=GENE_INDEX_ARRAY,
        variant_gene_array=VARIANT_GENE_ARRAY,
        variant_gene_schema=VARIANT_GENE_SCHEMA)

VARIANT_TRANSCRIPT_ARRAY = 'variant_transcript'
VARIANT_TRANSCRIPT_SCHEMA = """
  <notused: int8 not null>
  [chrom          = 1:25:0:1;
   pos            = 0:*:0:10000000;
   transcript_idx = 0:*:0:20]"""

VARIANT_TRANSCRIPT_STORE_QUERY = """
  store(
    redimension(
      index_lookup(
        apply(
          aio_input('{{path}}', 'num_attributes=3'),
          chrom,         int64(a0),
          pos,           int64(a1),
          transcript_id, a2,
          notused,       int8(0)) as INPUT,
        {transcript_index_array},
        INPUT.transcript_id,
        transcript_idx),
      {variant_transcript_schema}),
    {variant_transcript_array})""".format(
        transcript_index_array=TRANSCRIPT_INDEX_ARRAY,
        variant_transcript_array=VARIANT_TRANSCRIPT_ARRAY,
        variant_transcript_schema=VARIANT_TRANSCRIPT_SCHEMA)


# -- -
# -- - Load: COVERAGE - --
# -- -
COVERAGE_FILE = os.path.join(
    GBE_DATA_PATH, 'coverage', 'Panel2016.all.coverage.txt.gz')

COVERAGE_ARRAY = 'coverage'
COVERAGE_SCHEMA = """
  <odds_ratio:  double,
   log10pvalue: double,
   flag:        string,
   category:    string>
  [chrom = 1:25:0:1;
   pos   = 0:*:0:10000000]"""
COVERAGE_SCHEMA_OBJ = scidbpy.schema.Schema.fromstring(COVERAGE_SCHEMA)

COVERAGE_STORE_QUERY = """
  store(
    redimension(
      apply(
        filter(
          aio_input('{{path}}', 'num_attributes=8'),
          substr(a0, 0, 1) <> '#'),
        chrom,       int64(a0),
        pos,         int64(a1),
        odds_ratio,  dcast(a2, double(null)),
        log10pvalue, dcast(a5, double(null)),
        flag,        a6,
        category,    iif(a7 = 'transcript_ablation' or
                         a7 = 'splice_acceptor_variant' or
                         a7 = 'splice_donor_variant' or
                         a7 = 'stop_gained' or
                         a7 = 'frameshift_variant',
                         'lof_variant',
                         iif(a7 = 'stop_lost' or
                             a7 = 'start_lost' or
                             a7 = 'initiator_codon_variant' or
                             a7 = 'transcript_amplification' or
                             a7 = 'inframe_insertion' or
                             a7 = 'inframe_deletion' or
                             a7 = 'missense_variant',
                           'missense_variant',
                           iif(a7 = 'protein_altering_variant' or
                               a7 = 'splice_region_variant' or
                               a7 = 'incomplete_terminal_codon_variant' or
                               a7 = 'stop_retained_variant' or
                               a7 = 'synonymous_variant',
                             'synonymous_variant',
                             a7)))),
      {coverage_array_schema}),
    {coverage_array})""".format(coverage_array=COVERAGE_ARRAY,
                                coverage_array_schema=COVERAGE_SCHEMA)


# -- -
# -- - Load: DBSNP - --
# -- -
DBSNP_FILE = os.path.join(GBE_DATA_PATH, 'dbsnp150.txt.gz')

DBSNP_BY_RSID_ARRAY = 'dbsnp_by_rsid'
DBSNP_BY_RSID_SCHEMA = """
  <chrom: int64,
   pos:   int64>
  [rsid = 0:*:0:1000000]"""

DBSNP_BY_RSID_STORE_QUERY = """
  store(
    redimension(
      apply(
        aio_input('{{path}}', 'num_attributes=3'),
        rsid,  int64(a0),
        chrom, iif(a1 = 'X',
                   23,
                   iif(a1 = 'Y',
                       24,
                       dcast(a1, int64(0)))),
        pos,   int64(a2)),
      {dbsnp_by_rsid_schema}),
    {dbsnp_by_rsid_array})""".format(
        dbsnp_by_rsid_schema=DBSNP_BY_RSID_SCHEMA,
        dbsnp_by_rsid_array=DBSNP_BY_RSID_ARRAY)

DBSNP_BY_CHROM_POS_ARRAY = 'dbsnp_by_chrom_pos'
DBSNP_BY_CHROM_POS_SCHEMA = """
  <rsid: int64>
  [chrom = 1:25:0:1;
   pos   = 0:*:0:10000000]"""
DBSNP_BY_CHROM_POS_SCHEMA_OBJ = scidbpy.schema.Schema.fromstring(
    DBSNP_BY_CHROM_POS_SCHEMA)

DBSNP_BY_CHROM_POS_STORE_QUERY = """
  store(
    redimension({dbsnp_by_rsid_array},
                {dbsnp_by_chrom_pos_schema},
                false),
    {dbsnp_by_chrom_pos_array})""".format(
        dbsnp_by_rsid_array=DBSNP_BY_RSID_ARRAY,
        dbsnp_by_chrom_pos_schema=DBSNP_BY_CHROM_POS_SCHEMA,
        dbsnp_by_chrom_pos_array=DBSNP_BY_CHROM_POS_ARRAY)


# -- -
# -- - Load: BIM - --
# -- -
BIM_FILE = os.path.join(GBE_DATA_PATH, 'bims_combined.vep.cf.tsv.gz')

BIM_ARRAY = 'bim'
BIM_SCHEMA = """
  <ref:         string,
   alt:         string,
   consequence: string,
   hgvsp:       string>
  [chrom = 1:25:0:1;
   pos   = 0:*:0:10000000]"""

BIM_STORE_QUERY = """
  store(
    redimension(
      apply(
        aio_input('{{path}}', 'num_attributes=11', 'header=1'),
        chrom,       int64(a0),
        pos,         int64(a1),
        ref,         a2,
        alt,         a3,
        consequence, a8,
        hgvsp,       a9),
      {bim_schema}),
    {bim_array})""".format(
        bim_file=BIM_FILE,
        bim_schema=BIM_SCHEMA,
        bim_array=BIM_ARRAY)


# == =
# == = LOOKUP = ==
# == =
LOOKUP_QUERY = """
  cross_join(
    {main_array},
    filter({index_array}, {id_attr} = '{id_val}'),
    {main_array}.{idx_attr},
    {index_array}.{idx_attr})"""

EXISTS_QUERY = """
  aggregate(
    filter({array_name}, {attr_name} = {attr_val}),
    count(*))"""

EXISTS_SCHEMA = scidbpy.schema.Schema.fromstring('<count:int64>[notused]')

# -- -
# -- - Lookup: ICD - --
# -- -
ICD_INFO_MAP_QUERY = 'project({icd_info_array}, icd, Name)'.format(
    icd_info_array=ICD_INFO_ARRAY)

ICD_INFO_MAP_SCHEMA = scidbpy.schema.Schema.fromstring("""
  <icd:string, Name:string>[notused]""")

ICD_INFO_ICD_QUERY = 'project({icd_info_array}, icd)'.format(
    icd_info_array=ICD_INFO_ARRAY)

ICD_INFO_ICD_SCHEMA = scidbpy.schema.Schema.fromstring("""
  <icd:string>[notused]""")

ICD_INFO_IDX_MAX_QUERY = """
  aggregate(
    apply({icd_info_array}, idx, icd_idx),
     max(idx))""".format(icd_info_array=ICD_INFO_ARRAY)

ICD_INFO_IDX_MAX_SCHEMA = scidbpy.schema.Schema.fromstring("""
  <max:int64>[notused]""")

ICD_CHROM_POS_LOOKUP_QUERY = """
  equi_join(
    between({icd_array},
            null, {{chrom}}, {{start}}, null, null,
            null, {{chrom}}, {{stop}}, null, null),
    {{icd_info_filter}},
    'left_names=icd_idx',
    'right_names=icd_idx',
    'keep_dimensions=1',
    'algorithm=hash_replicate_right')""".format(icd_array=ICD_ARRAY)

ICD_CHROM_POS_LOOKUP_SCHEMA = scidbpy.schema.Schema.fromstring("""
  <icd_idx:     int64 not null,
   icdind:      int64,
   affyid:      string,
   or_val:      double,
   se:          double,
   pvalue:      double,
   lor:         double,
   log10pvalue: double,
   l95or:       double,
   u95or:       double,
   chrom:       int64 not null,
   pos:         int64 not null,
   pdecimal:    int64 not null,
   synthetic:   int64 not null,
   icd:         string,
   Case:        int64,
   Name:        string>
  [notused0;
   notused1]""")

ICD_AFFYID_LOOKUP_QUERY = """
  equi_join(
    {icd_affyid_array},
    filter({affyid_index_array}, affyid = '{{affyid}}'),
    'left_names=affyid_idx',
    'right_names=affyid_idx',
    'algorithm=hash_replicate_right')""".format(
        icd_affyid_array=ICD_AFFYID_ARRAY,
        affyid_index_array=AFFYID_INDEX_ARRAY)

ICD_AFFYID_LOOKUP_SCHEMA = scidbpy.schema.Schema.fromstring("""
  <affyid_idx: int64 not null,
   chrom:      int64,
   pos:        int64,
   affyid:     string>
  [notused0;
   notused1]""")

ICD_VARIANT_LOOKUP_QUERY = """
  equi_join(
    project({variant_array},
            rsid,
            ref,
            alt,
            filter,
            exac_nfe,
            csq),
    cross_join(
        project(
          between({icd_array},
                  null, null, null, {{pdecimal}}, null,
                  null, null, null, null, null),
          or_val,
          pvalue,
          log10pvalue),
        filter({icd_info_array}, icd = '{{icd}}'),
        {icd_array}.icd_idx,
        {icd_info_array}.icd_idx) as icd_join,
    'left_names=chrom,pos',
    'right_names=chrom,pos',
    'keep_dimensions=1',
    'algorithm=merge_right_first')""".format(
        icd_array=ICD_ARRAY,
        icd_info_array=ICD_INFO_ARRAY,
        variant_array=VARIANT_ARRAY)

ICD_VARIANT_SCAN_QUERY = """
  equi_join(
    project({variant_array},
            rsid,
            ref,
            alt,
            filter,
            exac_nfe,
            csq),
    cross_join(
        project(
          between({icd_array},
                  null, null, null, {{pdecimal}}, null,
                  null, null, null, null, null),
          or_val,
          pvalue,
          log10pvalue),
        filter({icd_info_array}, Case >= 100),
        {icd_array}.icd_idx,
        {icd_info_array}.icd_idx) as icd_join,
    'left_names=chrom,pos',
    'right_names=chrom,pos',
    'keep_dimensions=1',
    'algorithm=merge_right_first')""".format(
        icd_array=ICD_ARRAY,
        icd_info_array=ICD_INFO_ARRAY,
        variant_array=VARIANT_ARRAY)

VARIANT_X_ICD_X_INFO_SCHEMA = scidbpy.schema.Schema.fromstring("""
  <chrom:       int64 not null,
   pos:         int64 not null,
   rsid:        int64,
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
   Name:        string,
   icd_idx:     int64 not null,
   pdecimal:    int64 not null,
   synthetic:   int64 not null>
  [notused0;
   notused1]""")


# -- -
# -- - Lookup: GENE - --
# -- -
GENE_FILTER_QUERY = """
  equi_join(
    {{gene_array_filter}},
    {gene_index_array},
    'left_names=gene_idx',
    'right_names=gene_idx',
    'algorithm=hash_replicate_right')""".format(
        gene_index_array=GENE_INDEX_ARRAY)

GENE_FILTER_SCHEMA = scidbpy.schema.Schema.fromstring("""
  <gene_idx:       int64 not null,
   gene_name:      string,
   strand:         string,
   full_gene_name: string,
   omim_accession: string,
   transcript_idx: int64,
   chrom:          int64,
   start:          int64,
   stop:           int64,
   gene_id:        string>
  [notused0;
   notused1]""")

GENE_INDEX_LOOKUP_QUERY = """
  between({gene_index_array}, {{gene_idx}}, {{gene_idx}})""".format(
      gene_index_array=GENE_INDEX_ARRAY)

GENE_TRANSCRIPT_BY_ID_QUERY = """
  equi_join(
    equi_join({gene_array},
              filter({gene_index_array}, gene_id = '{{gene_id}}'),
              'left_names=gene_idx',
              'right_names=gene_idx',
              'algorithm=hash_replicate_right'),
    {transcript_index_array},
    'left_names=transcript_idx',
    'right_names=transcript_idx',
    'algorithm=hash_replicate_right')""".format(
        gene_array=GENE_ARRAY,
        gene_index_array=GENE_INDEX_ARRAY,
        transcript_index_array=TRANSCRIPT_INDEX_ARRAY)

GENE_TRANSCRIPT_BY_ID_SCHEMA = scidbpy.schema.Schema.fromstring("""
  <transcript_idx:    int64,
   gene_idx:          int64 not null,
   gene_name:         string,
   strand:            string,
   full_gene_name:    string,
   omim_accession:    string,
   chrom:             int64,
   start:             int64,
   stop:              int64,
   c_transcript_info: string,
   transcript_info:   string,
   exon_info:         string,
   gene_id:           string,
   transcript_id:     string>
  [notused0;
   notused1]""")

GENE_TRANSCRIPT_BY_IDX_QUERY = """
  equi_join(
    between({gene_array}, {{gene_idx}}, {{gene_idx}}),
    {transcript_index_array},
    'left_names=transcript_idx',
    'right_names=transcript_idx',
    'algorithm=hash_replicate_right')""".format(
        gene_array=GENE_ARRAY,
        transcript_index_array=TRANSCRIPT_INDEX_ARRAY)

GENE_TRANSCRIPT_BY_IDX_SCHEMA = scidbpy.schema.Schema.fromstring("""
  <transcript_idx:    int64,
   gene_name:         string,
   strand:            string,
   full_gene_name:    string,
   omim_accession:    string,
   chrom:             int64,
   start:             int64,
   stop:              int64,
   c_transcript_info: string,
   transcript_info:   string,
   exon_info:         string,
   transcript_id:     string>
  [notused0;
   notused1]""")

GENE_REGION_QUERY = """
  equi_join(
    filter({gene_array},
           chrom = {{chrom}} and start <= {{stop}} and stop >= {{start}}),
    {gene_index_array},
    'left_names=gene_idx',
    'right_names=gene_idx',
    'algorithm=hash_replicate_right')""".format(
        gene_array=GENE_ARRAY,
        gene_index_array=GENE_INDEX_ARRAY)

GENE_REGION_SCHEMA = scidbpy.schema.Schema.fromstring("""
  <gene_idx:          int64 not null,
   gene_name:         string,
   strand:            string,
   full_gene_name:    string,
   omim_accession:    string,
   transcript_idx:    int64,
   chrom:             int64,
   start:             int64,
   stop:              int64,
   c_transcript_info: string,
   transcript_info:   string,
   exon_info:         string,
   gene_id:           string>
  [notused0;
   notused1]""")

GENE_ID_BY_NAME_QUERY = """
  project(
    equi_join(
      filter({gene_array}, gene_name = '{{gene_name}}'),
      {gene_index_array},
      'left_names=gene_idx',
      'right_names=gene_idx',
      'algorithm=hash_replicate_right'),
    gene_id)""".format(gene_array=GENE_ARRAY,
                       gene_index_array=GENE_INDEX_ARRAY)

GENE_ID_BY_NAME_SCHEMA = scidbpy.schema.Schema.fromstring(
    '<gene_id:string>[notused]')

GENE_VARIANT_LOOKUP = """
  sort(
    equi_join(
      equi_join(
        project({variant_array}, rsid, ref, alt, exac_nfe),
        project(
          equi_join(
            {variant_gene_array},
            equi_join(
              {gene_index_array},
              project({{gene_filter}}, gene_name),
              'left_names=gene_idx',
              'right_names=gene_idx',
              'algorithm=hash_replicate_right'),
            'left_names=gene_idx',
            'right_names=gene_idx',
            'keep_dimensions=1',
            'algorithm=hash_replicate_right'),
          chrom,
          pos,
          gene_name),
        'left_names=chrom,pos',
        'right_names=chrom,pos'),
      {bim_array},
      'left_names=chrom,pos,ref,alt',
      'right_names=chrom,pos,ref,alt',
      'algorithm=hash_replicate_right'),
    chrom,
    pos,
    ref,
    alt)""".format(variant_array=VARIANT_ARRAY,
                   variant_gene_array=VARIANT_GENE_ARRAY,
                   gene_index_array=GENE_INDEX_ARRAY,
                   bim_array=BIM_ARRAY)

GENE_VARIANT_SCHEMA = scidbpy.schema.Schema.fromstring("""
  <chrom:       int64 not null,
   pos:         int64 not null,
   ref:         string,
   alt:         string,
   rsid:        int64,
   exac_nfe:    double,
   gene_name:   string,
   consequence: string,
   hgvsp:       string>
  [notused0;
   notused01]""")

VARIANT_ICD_LOOKUP = """
  project(
    sort(
      equi_join(
        project(
          equi_join(
            project({icd_array}, se, pvalue, lor),
            project({{icd_filter}}, Case),
            'left_names=icd_idx',
            'right_names=icd_idx',
            'keep_dimensions=1',
            'algorithm=hash_replicate_right'),
          se,
          lor,
          chrom,
          pos),
        project(
          equi_join(
            project({variant_array}, ref, alt, rsid),
            project(
              equi_join(
                {variant_gene_array},
                project(
                  equi_join(
                    {gene_index_array},
                    project({{gene_filter}}, chrom),
                    'left_names=gene_idx',
                    'right_names=gene_idx',
                    'algorithm=hash_replicate_right'),
                  gene_idx),
                'left_names=gene_idx',
                'right_names=gene_idx',
                'keep_dimensions=1',
                'algorithm=hash_replicate_right'),
              chrom,
              pos,
              gene_idx),
            'left_names=chrom,pos',
            'right_names=chrom,pos'),
          chrom,
          pos, ref, alt),
        'left_names=chrom,pos',
        'right_names=chrom,pos',
        'right_outer=1'),
      chrom,
      pos,
      ref,
      alt),
    se,
    lor)""".format(icd_array=ICD_ARRAY,
                   variant_array=VARIANT_ARRAY,
                   variant_gene_array=VARIANT_GENE_ARRAY,
                   gene_index_array=GENE_INDEX_ARRAY)

VARIANT_ICD_SCHEMA = scidbpy.schema.Schema.fromstring("""
  <se:     double,
   lor:    double>
  [notused0;
   notused1]""")

TRANSCRIPT_OR_EXON_BETWEEN_QUERY = """
  between({array_name},
          {gene_idx}, {transcript_idx}, null, null, null, null,
          {gene_idx}, {transcript_idx}, null, null, null, null)"""

TRANSCRIPT_ID_LOOKUP = """
  cross_join(
    between({transcript_array},
            {{gene_idx}}, null,
            {{gene_idx}}, null),
    {transcript_index_array},
    {transcript_array}.transcript_idx,
    {transcript_index_array}.transcript_idx)""".format(
        transcript_array=TRANSCRIPT_ARRAY,
        transcript_index_array=TRANSCRIPT_INDEX_ARRAY)

TRANSCRIPT_ID_SCHEMA = scidbpy.schema.Schema.fromstring(
    TRANSCRIPT_SCHEMA.replace('>', ',{}>'.format(
            TRANSCRIPT_INDEX_SCHEMA[TRANSCRIPT_INDEX_SCHEMA.index('<') + 1:
                                    TRANSCRIPT_INDEX_SCHEMA.index('>')])))

TRANSCRIPT_GENE_LOOKUP = """
  cross_join(
    cross_join(
      {transcript_array},
      filter({transcript_index_array}, transcript_id = '{{transcript_id}}'),
      {transcript_array}.transcript_idx,
      {transcript_index_array}.transcript_idx),
    {gene_index_array},
    {transcript_array}.gene_idx,
    {gene_index_array}.gene_idx)""".format(
        transcript_array=TRANSCRIPT_ARRAY,
        transcript_index_array=TRANSCRIPT_INDEX_ARRAY,
        gene_index_array=GENE_INDEX_ARRAY)

TRANSCRIPT_GENE_SCHEMA = scidbpy.schema.Schema.fromstring(
    TRANSCRIPT_SCHEMA.replace(
        '>',
        ',{}, {}>'.format(
            TRANSCRIPT_INDEX_SCHEMA[TRANSCRIPT_INDEX_SCHEMA.index('<') + 1:
                                    TRANSCRIPT_INDEX_SCHEMA.index('>')],
            GENE_INDEX_SCHEMA[GENE_INDEX_SCHEMA.index('<') + 1:
                              GENE_INDEX_SCHEMA.index('>')])))


# -- -
# -- - Lookup: VARIANT - --
# -- -
VARIANT_LOOKUP_QUERY = """
  between({variant_array},
          {{chrom}}, {{start}},
          {{chrom}}, {{stop}})""".format(
    variant_array=VARIANT_ARRAY)

VARIANT_LOOKUP_SCHEMA = scidbpy.schema.Schema.fromstring(VARIANT_SCHEMA)

VARIANT_LIMIT_QUERY = """
  limit(
    project(
      between({variant_array},
              {{chrom}}, {{start}},
              {{chrom}}, {{start}}),
      rsid,
      ref,
      alt,
      filter,
      csq),
    1)""".format(
    variant_array=VARIANT_ARRAY)

VARIANT_LIMIT_SCHEMA = scidbpy.schema.Schema.fromstring("""
  <rsid:         int64,
   ref:          string,
   alt:          string,
   filter:       string,
   csq:          string>
  [chrom = 1:25:0:1;
   pos   = 0:*:0:10000000]""")

VARIANT_CHROM_POS_BY_RSID_QUERY = """
  limit(
    project(
      filter({variant_array}, rsid = {{rsid}}),
      rsid),
    2)""".format(variant_array=VARIANT_ARRAY)

VARIANT_CHROM_POS_BY_RSID_SCHEMA = scidbpy.schema.Schema.fromstring("""
  <rsid:int64>
  [chrom = 1:25:0:1;
   pos   = 0:*:0:10000000]""")

VARIANT_GENE_LOOKUP = """
  cross_join(
    {variant_array},
    between({variant_gene_array},
            null, null, {{gene_idx}},
            null, null, {{gene_idx}}),
    {variant_array}.chrom,
    {variant_gene_array}.chrom,
    {variant_array}.pos,
    {variant_gene_array}.pos)""".format(
        variant_array=VARIANT_ARRAY,
        variant_gene_array=VARIANT_GENE_ARRAY)

VARIANT_GENE_SCHEMA = scidbpy.schema.Schema.fromstring(
    VARIANT_GENE_SCHEMA.replace(
        '<',
        '<{},'.format(
            VARIANT_SCHEMA[VARIANT_SCHEMA.index('<') + 1:
                           VARIANT_SCHEMA.index('>')])))

VARIANT_TRANSCRIPT_IDX_LOOKUP = """
  cross_join(
    {variant_array},
    between({variant_transcript_array},
            null, null, {{transcript_idx}},
            null, null, {{transcript_idx}}),
    {variant_array}.chrom,
    variant_transcript.chrom,
    {variant_array}.pos,
    variant_transcript.pos)""".format(
        variant_array=VARIANT_ARRAY,
        variant_transcript_array=VARIANT_TRANSCRIPT_ARRAY)

VARIANT_X_TRANSCRIPT_SCHEMA = scidbpy.schema.Schema.fromstring(
    VARIANT_TRANSCRIPT_SCHEMA.replace(
        '<',
        '<{},'.format(
            VARIANT_SCHEMA[VARIANT_SCHEMA.index('<') + 1:
                           VARIANT_SCHEMA.index('>')])))

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

# -- -
# -- - Lookup: COVERAGE - --
# -- -
COVERAGE_LOOKUP_QUERY = """
  between({coverage_array},
          {{chrom_start}}, {{pos_start}},
          {{chrom_stop}}, {{pos_stop}})""".format(
                                coverage_array=COVERAGE_ARRAY)


# -- -
# -- - Lookup: DBSNP - --
# -- -
DBSNP_LOOKUP_QUERY = """
  between({dbsnp_by_chrom_pos_array},
          {{chrom}}, {{pos}},
          {{chrom}}, {{pos}})""".format(
              dbsnp_by_chrom_pos_array=DBSNP_BY_CHROM_POS_ARRAY)

DBSNP_VARIANT_LOOKUP_QUERY = """
  equi_join(
    between({dbsnp_by_rsid_array}, {{rsid}}, {{rsid}}),
    project({variant_array}, rsid),
    'left_names=chrom,pos',
    'right_names=chrom,pos',
    'algorithm=hash_replicate_left')""".format(
        dbsnp_by_rsid_array=DBSNP_BY_RSID_ARRAY,
        variant_array=VARIANT_ARRAY)

DBSNP_VARIANT_LOOKUP_SCHEMA = scidbpy.schema.Schema.fromstring("""
  <chrom:  int64,
   pos:    int64,
   rsid:   int64>
  [notused0;
   notused1]""")
