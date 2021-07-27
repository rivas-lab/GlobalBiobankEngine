import itertools
import re
#import scidbpy
import scidbbiobank 
import config
import utils
import logging
import numpy 
import pandas 
logging.basicConfig(filename='loglookups.txt', level=logging.DEBUG)
UNSUPPORTED_QUERIES = set((
    'CMD1G'
#    'CMH9',
#    'CMPD4',
#    'ENSG00000155657',
#    'ENST00000342175',
#    'ENST00000342992',
#    'ENST00000359218',
#    'ENST00000460472',
#    'ENST00000589042',
#    'ENST00000591111',
#    'FLJ32040',
#    'LGMD2J',
#    'MYLK5',
#    'TMD',
#    'TTN',
))

RSID_FORMAT = '{chrom}-{pos}-{ref}-{alt}'

# 1:1-1000
REGION_RE1 = re.compile(r'^(\d+|X|Y|M|MT)\s*:\s*(\d+)-(\d+)$')
REGION_RE2 = re.compile(r'^(\d+|X|Y|M|MT)\s*:\s*(\d+)$')
REGION_RE3 = re.compile(r'^(\d+|X|Y|M|MT)$')
REGION_RE4 = re.compile(r'^(\d+|X|Y|M|MT)\s*[-:]\s*(\d+)-([ATCG]+)-([ATCG]+)$')

TRANSCRIPT_INFO_KEYS = ('transcript_id', 'strand', 'chrom', 'start', 'stop')
EXON_INFO_KEYS = ('feature_type', 'chrom', 'start', 'stop')


def numpy2dict0(ar):
    """Convert SciDB NumPy record result to Python dictionary and populate
    nullable attributes with values (discards null codes).

    """
    if not len(ar):
        return None
    el = ar[0]
    return dict(
        (de[0],
         el[de[0]]['val'] if isinstance(de[1], list) else el[de[0]])
        for de in ar.dtype.descr if de[0] != 'notused')


def numpy2dict(ar):
    """Convert SciDB NumPy array result to a list of Python dictionary and
    populate nullable attributes with values (discards null codes).

    """
    return [
        dict(
            (de[0],
             el[de[0]]['val'] if isinstance(de[1], list) else el[de[0]])
            for de in ar.dtype.descr if de[0] != 'notused'
        )
        for el in ar]




def parse_vep_annotations(csq, gene_id=None, transcript_id=None):
    return [ann for ann in (dict(zip(config.VARIANT_CSQ, cs.split('|')))
                            for cs in csq.split(','))
            if ('Feature' in ann and
                ann['Feature'].startswith('ENST') and
                (gene_id is None or ann['Gene'] == gene_id) and
                (transcript_id is None or ann['Feature'] == transcript_id))]


def format_variants(variants, add_ann=False, gene_id=None, transcript_id=None):
    for variant in variants:
        variant['rsid'] = ('rs{}'.format(variant['rsid'])
                           if variant['rsid'] else '.')
        variant['variant_id'] = RSID_FORMAT.format(
            chrom=variant['chrom'],
            pos=variant['pos'],
            ref=variant['ref'],
            alt=variant['alt'])

        vep_annotations = parse_vep_annotations(
            variant['csq'], gene_id, transcript_id)
        if add_ann:
            variant['vep_annotations'] = vep_annotations

        variant['genes'] = list(set(ann['Gene'] for ann in vep_annotations))
        variant['gene_name'] = ','.join(variant['genes'][:3])
        variant['gene_symbol'] = ','.join(
            itertools.islice(set(ann['SYMBOL'] for ann in vep_annotations), 3))
        variant['transcripts'] = list(set(
            ann['Feature'] for ann in vep_annotations))

        utils.add_consequence_to_variant(variant, vep_annotations)

    return variants


def cast_pos_info(gene):
    for key in ('chrom', 'start', 'stop'):
        if key in gene:
            gene[key] = int(gene[key])
    return gene


def add_xpos(gene):
    if gene and all(k in gene.keys() for k in ('chrom', 'start', 'end')):
        gene['xstart'] = gene['chrom'] * config.XOFF + gene['start']
        gene['xstop'] = gene['chrom'] * config.XOFF + gene['end']
    return gene


def exists(db, array_name, attr_name, attr_val):
    """
    Search bar

    MongoDB:
      db.genes.find({'gene_id': 'ENSG00000107404'}, fields={'_id': False})

    SciDB:
      aggregate(
        filter(gene_index, gene_id = 'ENSG00000198734'),
        count(*));
    """
    return bool(
        db.iquery(
            config.EXISTS_QUERY.format(
                array_name=array_name,
                attr_name=attr_name,
                attr_val=attr_val),
            schema=config.EXISTS_SCHEMA,
            fetch=True,
            atts_only=True)[0]['count']['val'])


# -- -
# -- - ICD - --
# -- -
def get_icd_name_map(db):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/coding/RH117

    MongoDB:
      db.icd_info.find({'icd': 'RH117'}, fields={'_id': False})

    SciDB:
      project(icd_info, icd, Name);
    """
    return dict((i['icd']['val'], '&nbsp;'.join(i['Name']['val'].split()))
                for i in db.iquery(config.ICD_INFO_MAP_QUERY,
                                   schema=config.ICD_INFO_MAP_SCHEMA,
                                   fetch=True,
                                   atts_only=True))


def exists_icd(db, icd):
    """
    Search bar

    MongoDB:
      db.icd_info.find({'icd': 'RH141'}, fields={'_id': False})

    SciDB:
      aggregate(
        filter(icd_info, icd = 'RH141'),
        count(*));

    SciDBnew: 
    res = db.get_phenotype_fields(association_set=str(db.list_association_sets()['name'][0]))
    resphe = [res['description'] == icd]
    return bool(resphe.empty)
    """
    res = db.get_phenotype_fields(association_set=str(db.list_association_sets()['name'][0]))
    resphe = res[res['description'] == icd]
    return not bool(resphe.empty)

def get_phe_title(db, phename):
    """
    Search bar
    SciDBnew: 
    res = str(list(phef[phef['description'] == "asthma_diagnosed_by_doctor"]['title'])[0])
    """
    phef = db.get_phenotype_fields(association_set=str(db.list_association_sets()['name'][0]))
    if phef.empty:
        return None
    else:
        res = str(list(phef[phef['description'] == phename]['title'])[0])
    return res


def get_phe_name(db, icd):
    """
    Search bar
    SciDBnew: 

    icdres['shortname'] = icdres.apply(lambda row: str(df[df['title'] == row['title']]['notes'].squeeze()).split(';')[1].split('=')[1], axis = 1)
    logging.info('shortname')
    icdres['Case'] = icdres.apply(lambda row: str(df[df['title'] == row['title']]['notes'].squeeze()).split(';')[0].split('=')[1], axis = 1)
    icdres['gene_name'], icdres['gene_symbol'], icdres['HGVSp'], icdres['HGVSc'] = zip(*icdres['annotations'].map(utils.return_gene_vep))
    #icdres.apply(lambda row: utils.return_gene_vep(row['annotations']), axis = 1)                                                                                                
    logging.info(icdres['gene_symbol'])
    logging.info('Case')
    icdres['log10pvalue'] = icdres.apply(lambda row: -numpy.log10(row['pvalue']), axis = 1)
    logging.info('l10pval')
    icdres['icd'] = icdres['title']
    logging.info('icd')
    if icdres.loc[icdres['odds_ratio'].isna()].shape[0] < icdres.loc[icdres['beta'].isna()].shape[0]:
        icdres['or_val'] = icdres['odds_ratio']

    res = str(list(phef[phef['description'] == "asthma_diagnosed_by_doctor"]['title'])[0])
    """
    phef = db.get_phenotype_fields(association_set=str(db.list_association_sets()['name'][0]))
    if phef.empty:
        return None
    else:
        res = str(phef[phef['title'] == icd]['notes'].squeeze()).split(';')[1].split('=')[1]
    return res


def get_phe_case(db, icd):
    """
    Search bar
    SciDBnew: 

    icdres['shortname'] = icdres.apply(lambda row: str(df[df['title'] == row['title']]['notes'].squeeze()).split(';')[1].split('=')[1], axis = 1)
    logging.info('shortname')
    icdres['Case'] = icdres.apply(lambda row: str(df[df['title'] == row['title']]['notes'].squeeze()).split(';')[0].split('=')[1], axis = 1)
    icdres['gene_name'], icdres['gene_symbol'], icdres['HGVSp'], icdres['HGVSc'] = zip(*icdres['annotations'].map(utils.return_gene_vep))
    #icdres.apply(lambda row: utils.return_gene_vep(row['annotations']), axis = 1)                                                                                                
    logging.info(icdres['gene_symbol'])
    logging.info('Case')
    icdres['log10pvalue'] = icdres.apply(lambda row: -numpy.log10(row['pvalue']), axis = 1)
    logging.info('l10pval')
    icdres['icd'] = icdres['title']
    logging.info('icd')
    if icdres.loc[icdres['odds_ratio'].isna()].shape[0] < icdres.loc[icdres['beta'].isna()].shape[0]:
        icdres['or_val'] = icdres['odds_ratio']

    res = str(list(phef[phef['description'] == "asthma_diagnosed_by_doctor"]['title'])[0])
    """
    phef = db.get_phenotype_fields(association_set=str(db.list_association_sets()['name'][0]))
    if phef.empty:
        return None
    else:
        res = str(phef[phef['title'] == icd]['notes'].squeeze()).split(';')[0].split('=')[1]
    return res


def get_icd_by_chrom_pos(db, chrom, start, stop=None, icd=None):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/variant/1-39381448

    MongoDB:
      db.icd.find({'xpos': '1039381448'}, fields={'_id': False})

    SciDB:
      equi_join(
        between(icd, null, 1, 39381448, null,
                     null, 1, 39381448, null),
        icd_info,
        'left_names=icd_idx',
        'right_names=icd_idx',
        'keep_dimensions=1',
        'algorithm=hash_replicate_right');
    SciDBnew:
    db.get_association_data(association_set=assocset, chromosome=22, position = 32334104)
    """
    if not stop:
        stop = start
#    if not icd:
#        icd_info_filter = config.ICD_INFO_ARRAY
#    else:
#        icd_info_filter = 'filter({icd_info_array}, {cond})'.format(
#            icd_info_array=config.ICD_INFO_ARRAY,
#            cond=' or '.join("icd = '{}'".format(i) for i in icd))
#    return numpy2dict(
#        db.iquery(
#            config.ICD_CHROM_POS_LOOKUP_QUERY.format(
#                icd_info_filter=icd_info_filter,
#                chrom=chrom,
#                start=start,
#                stop=stop),
#            schema=config.ICD_CHROM_POS_LOOKUP_SCHEMA,
#            fetch=True,
#            atts_only=True))


def get_icd_affyid(db, affyid):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/intensity/723307

    MongoDB:
      db.icd.find_one({'affyid': '723307'}, fields={'_id': False})

    SciDB:
      equi_join(
        icd_affyid,
        filter(affyid_index, affyid = '723307'),
        'left_names=affyid_idx',
        'right_names=affyid_idx',
        'algorithm=hash_replicate_right');
    """
    return numpy2dict0(
        db.iquery(
            config.ICD_AFFYID_LOOKUP_QUERY.format(affyid=affyid),
            schema=config.ICD_AFFYID_LOOKUP_SCHEMA,
            fetch=True,
            atts_only=True))


def get_icd_variant_by_icd_id_pvalue(db, icd_id, field_identifier, pvalue=0.001):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/coding/RH117

    MongoDB:
      db.icd.find({'icd': 'RH117', 'stats.pvalue': {'$lt': 0.01}},
                  fields={'_id': false})
      db.icd_info.find({'icd': 'RH117'}, fields={'_id': False})
      db.variants.find({'xpos': '1039381448'}, fields={'_id': False})

    SciDB:
      equi_join(
        project(variant,
                rsid,
                ref,
                alt,
                filter,
                exac_nfe,
                csq),
        cross_join(
            project(
              between(icd, null, null, null, 1,    null,
                           null, null, null, null, null),
              or_val,
              pvalue,
              log10pvalue),
            filter(icd_info, icd = 'RH117'),
            icd.icd_idx,
            icd_info.icd_idx) as icd_join,
        'left_names=chrom,pos',
        'right_names=chrom,pos',
        'keep_dimensions=1',
        'algorithm=merge_right_first');

    SciDBnew: 
    bb.get_association_data(association_set=assocset,
          chromosome=22, start=32300000, end=32400000, pvalue_max=1, field_id = int(df[df['title'] == 'HC382']['field_id']))   
    """
    assocset = str(db.list_association_sets()['name'][0])
    df = db.get_phenotype_fields(association_set=assocset, include_pvalue_threshold=True)
    if 'enrichlogor' in list(db.get_variant_fields()['name']):
        icdres = db.get_association_data(association_set=assocset, pvalue_max = pvalue, field_id = field_identifier, variant_fields = ('chrom','pos','ref','alt', 'variant_identity', 'major_consequence', 'category', 'gene_name', 'gene_symbol', 'HGVSp', 'HGVSc',  'consequence', 'all_filters', 'maf', 'ukbb_freq', 'ld', 'rsid', 'enrichp', 'enrichlogor'), association_fields = ('pvalue', 'title', 'beta', 'odds_ratio', 'beta', 'se'), se_range=(0, 0.21))
    else:
        icdres = db.get_association_data(association_set=assocset, pvalue_max = pvalue, field_id = field_identifier, variant_fields = ('chrom','pos','ref','alt', 'variant_identity', 'major_consequence', 'category', 'gene_name', 'gene_symbol', 'HGVSp', 'HGVSc',  'consequence', 'all_filters',  'maf', 'ukbb_freq', 'ld', 'rsid'), association_fields = ('pvalue', 'title', 'beta', 'odds_ratio', 'beta', 'se'), se_range=(0, 0.21))
   # icdres = db.get_association_data(association_set=assocset, pvalue_max = pvalue, variant_fieldsfield_id = field_identifier, association_fields = ('pvalue', 'title', 'beta', 'odds_ratio', 'beta', 'se'))
#    icdres = db.get_association_data(association_set=assocset, field_id = field_identifier)
    ## SE filter
#    icdres['variant_id'] = icdres.apply(lambda row: RSID_FORMAT.format(
#            chrom=row['chrom'],
#            pos=row['pos'],
#            ref=row['ref'],
#            alt=row['alt']), axis = 1)
#    icdres.loc[~icdres.consequence.isin(utils.csq_order),'consequence'] = 'intergenic'
#    icdres['major_consequence'] = icdres['consequence']
#    icdres['category'] = icdres.apply(lambda row: utils.add_category_to_variant(row['consequence']), axis = 1)
#    icdres = icdres[icdres.se <= .2] 
    icdres['filter'] = icdres.apply(lambda row: 'PASS' if (row['all_filters'] != 1 and row['all_filters'] != 2 and row['all_filters'] != 3) else 'FAIL', axis = 1)
    icdres['pvalue'] = icdres.apply(lambda row: row['pvalue'] if row['pvalue'] != 0 else 1e-300, axis = 1)
   # icdres['Name'] =  icdres.apply(lambda row: df[df['title'] == row['title']]['description'].squeeze(), axis = 1)
    icdres['shortname'] = icdres.apply(lambda row: str(df[df['title'] == row['title']]['notes'].squeeze()).split(';')[1].split('=')[1], axis = 1)
    icdres['Name'] = icdres['shortname']
    icdres['Case'] = icdres.apply(lambda row: str(df[df['title'] == row['title']]['notes'].squeeze()).split(';')[0].split('=')[1], axis = 1)
#    icdres['gene_name'], icdres['gene_symbol'], icdres['HGVSp'], icdres['HGVSc'] = zip(*icdres['annotations'].map(utils.return_gene_vep)) 
   # icdres = icdres.drop(columns=['annotations'])
    #icdres.apply(lambda row: utils.return_gene_vep(row['annotations']), axis = 1)
    icdres['log10pvalue'] = icdres.apply(lambda row: -numpy.log10(row['pvalue']), axis = 1)
    icdres['icd'] = icdres['title']
    if icdres.loc[icdres['odds_ratio'].isna()].shape[0] < icdres.loc[icdres['beta'].isna()].shape[0]:
        icdres['or_val'] = icdres['odds_ratio']
        icdres['lor_val'] = icdres.apply(lambda row: numpy.log(row['odds_ratio']), axis = 1)
    else:
        icdres['or_val'] = icdres['beta']
        icdres['lor_val'] = icdres['beta']
 #   icdres['ukbb_freq'] = icdres['maf']
  #  logging.info(icdres)
    return icdres





def get_icd_variant_by_beta_pvalue(db, pvalue=0.00000005, betaabs = -.2):
    """
    """
    assocset = str(db.list_association_sets()['name'][0])
    df = db.get_phenotype_fields(association_set=assocset, include_pvalue_threshold=True)
    if 'enrichlogor' in list(db.get_variant_fields()['name']):
        icdres = db.get_association_data(association_set=assocset, pvalue_max = pvalue, variant_fields = ('chrom','pos','ref','alt', 'variant_identity', 'major_consequence', 'category', 'gene_name', 'gene_symbol', 'HGVSp', 'HGVSc',  'consequence', 'all_filters', 'maf', 'ukbb_freq', 'ld', 'rsid', 'enrichp', 'enrichlogor'), association_fields = ('pvalue', 'title', 'beta', 'odds_ratio', 'beta', 'se'), se_range=(0, 0.21), beta_range = (None, betaabs))
    else:
        icdres = db.get_association_data(association_set=assocset, pvalue_max = pvalue, variant_fields = ('chrom','pos','ref','alt', 'variant_identity', 'major_consequence', 'category', 'gene_name', 'gene_symbol', 'HGVSp', 'HGVSc',  'consequence', 'all_filters',  'maf', 'ukbb_freq', 'ld', 'rsid'), association_fields = ('pvalue', 'title', 'beta', 'odds_ratio', 'beta', 'se'), se_range=(0, 0.21), beta_range = (None, betaabs))
    if 'enrichlogor' in list(db.get_variant_fields()['name']):
        icdres2 = db.get_association_data(association_set=assocset, pvalue_max = pvalue, variant_fields = ('chrom','pos','ref','alt', 'variant_identity', 'major_consequence', 'category', 'gene_name', 'gene_symbol', 'HGVSp', 'HGVSc',  'consequence', 'all_filters', 'maf', 'ukbb_freq', 'ld', 'rsid', 'enrichp', 'enrichlogor'), association_fields = ('pvalue', 'title', 'beta', 'odds_ratio', 'beta', 'se'), se_range=(0, 0.21), beta_range = (-(betaabs), None))
    else:
        icdres2 = db.get_association_data(association_set=assocset, pvalue_max = pvalue, variant_fields = ('chrom','pos','ref','alt', 'variant_identity', 'major_consequence', 'category', 'gene_name', 'gene_symbol', 'HGVSp', 'HGVSc',  'consequence', 'all_filters',  'maf', 'ukbb_freq', 'ld', 'rsid'), association_fields = ('pvalue', 'title', 'beta', 'odds_ratio', 'beta', 'se'), se_range=(0, 0.21), beta_range = (-betaabs, None))
    icdres = icdres.append(icdres2)
   # icdres = db.get_association_data(association_set=assocset, pvalue_max = pvalue, variant_fieldsfield_id = field_identifier, association_fields = ('pvalue', 'title', 'beta', 'odds_ratio', 'beta', 'se'))
#    icdres = db.get_association_data(association_set=assocset, field_id = field_identifier)
    ## SE filter
#    icdres['variant_id'] = icdres.apply(lambda row: RSID_FORMAT.format(
#            chrom=row['chrom'],
#            pos=row['pos'],
#            ref=row['ref'],
#            alt=row['alt']), axis = 1)
#    icdres.loc[~icdres.consequence.isin(utils.csq_order),'consequence'] = 'intergenic'
#    icdres['major_consequence'] = icdres['consequence']
#    icdres['category'] = icdres.apply(lambda row: utils.add_category_to_variant(row['consequence']), axis = 1)
#    icdres = icdres[icdres.se <= .2] 
    icdres['filter'] = icdres.apply(lambda row: 'PASS' if (row['all_filters'] != 1 and row['all_filters'] != 2 and row['all_filters'] != 3) else 'FAIL', axis = 1)
    icdres['pvalue'] = icdres.apply(lambda row: row['pvalue'] if row['pvalue'] != 0 else 1e-300, axis = 1)
   # icdres['Name'] =  icdres.apply(lambda row: df[df['title'] == row['title']]['description'].squeeze(), axis = 1)
    icdres['shortname'] = icdres.apply(lambda row: str(df[df['title'] == row['title']]['notes'].squeeze()).split(';')[1].split('=')[1], axis = 1)
    icdres['Name'] = icdres['shortname']
    cond = icdres['shortname'].str.contains('Diet') | icdres['shortname'].str.contains('drinking habits') | icdres['shortname'].str.contains('bone mineral density') | icdres['shortname'].str.contains('drinking session') 
    icdres = icdres.drop(icdres[cond].index.values)
    icdres['Case'] = icdres.apply(lambda row: str(df[df['title'] == row['title']]['notes'].squeeze()).split(';')[0].split('=')[1], axis = 1)
#    icdres['gene_name'], icdres['gene_symbol'], icdres['HGVSp'], icdres['HGVSc'] = zip(*icdres['annotations'].map(utils.return_gene_vep)) 
   # icdres = icdres.drop(columns=['annotations'])
    #icdres.apply(lambda row: utils.return_gene_vep(row['annotations']), axis = 1)
    icdres['log10pvalue'] = icdres.apply(lambda row: -numpy.log10(row['pvalue']), axis = 1)
    icdres['icd'] = icdres['title']
    if icdres.loc[icdres['odds_ratio'].isna()].shape[0] < icdres.loc[icdres['beta'].isna()].shape[0]:
        icdres['or_val'] = icdres['odds_ratio']
        icdres['lor_val'] = icdres.apply(lambda row: numpy.log(row['odds_ratio']), axis = 1)
    else:
        icdres['or_val'] = icdres['beta']
        icdres['lor_val'] = icdres['beta']
 #   icdres['ukbb_freq'] = icdres['maf']
  #  logging.info(icdres)
    return icdres

def get_icd_variant_by_icd_id_pvalue_uncertain(db, icd_id, pvalue=0.001):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/coding/RH117

    MongoDB:
      db.icd.find({'icd': 'RH117', 'stats.pvalue': {'$lt': 0.01}},
                  fields={'_id': false})
      db.icd_info.find({'icd': 'RH117'}, fields={'_id': False})
      db.variants.find({'xpos': '1039381448'}, fields={'_id': False})

    SciDB:
      equi_join(
        project(variant,
                rsid,
                ref,
                alt,
                filter,
                exac_nfe,
                csq),
        cross_join(
            project(
              between(icd, null, null, null, 1,    null,
                           null, null, null, null, null),
              or_val,
              pvalue,
              log10pvalue),
            filter(icd_info, icd = 'RH117'),
            icd.icd_idx,
            icd_info.icd_idx) as icd_join,
        'left_names=chrom,pos',
        'right_names=chrom,pos',
        'keep_dimensions=1',
        'algorithm=merge_right_first');
    """
    pdecimal = config.ICD_PVALUE_MAP.get(pvalue, 0)
    return format_variants(
        numpy2dict(
            db.iquery(
                config.ICD_VARIANT_LOOKUP_UNCERTAIN_QUERY.format(
                    icd=icd_id, pdecimal=pdecimal),
                schema=config.VARIANT_X_ICD_X_INFO_SCHEMA,
                fetch=True,
                atts_only=True)))


def get_icd_variant_by_pvalue(db, pvalue=0.001,orvalue=0):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/target/1

    MongoDB:
      db.icd.find({'stats.pvalue': {'$lt': 0.01}}, fields={'_id': false})
      db.variants.find({'xpos': '1039381448'}, fields={'_id': False})

    SciDB:
      equi_join(
        project(variant,
                rsid,
                ref,
                alt,
                filter,
                exac_nfe,
                csq),
        cross_join(
            project(
              between(icd, null, null, null, 1,    null,
                           null, null, null, null, null),
              or_val,
              pvalue,
              log10pvalue),
            filter(icd_info, Cases >= 500),
            icd.icd_idx,
            icd_info.icd_idx) as icd_join,
        'left_names=chrom,pos',
        'right_names=chrom,pos',
        'keep_dimensions=1',
        'algorithm=merge_right_first');
    """
    pdecimal = config.ICD_PVALUE_MAP.get(pvalue, 0)
    print(pvalue,pdecimal)
    return format_variants(
        numpy2dict(
            db.iquery(
                config.ICD_VARIANT_SCAN_QUERY.format(pdecimal=pdecimal,orvalue=orvalue,pvalue=pvalue),
                schema=config.VARIANT_X_ICD_X_INFO_SCHEMA,
                fetch=True,
                atts_only=True)))


# -- -
# -- - GENE - --
# -- -
def get_gene_by_gene_names(db, gene_names=None):
    """
    SciDB:
      equi_join(
        gene,
        gene_index,
        'left_names=gene_idx',
        'right_names=gene_idx',
        'algorithm=hash_replicate_right');

      equi_join(
        filter(gene, gene_name = 'RP11-126K1.9' or gene_name='RP4-621B10.8'),
        gene_index,
        'left_names=gene_idx',
        'right_names=gene_idx',
        'algorithm=hash_replicate_right');
    """
    if gene_names:
        query = config.GENE_FILTER_QUERY.format(
            gene_array_filter='filter({gene_array}, {cond})'.format(
                gene_array=config.GENE_ARRAY,
                cond=' or '.join(
                    "gene_name = '{}'".format(g) for g in gene_names)))
    else:
        query = config.GENE_FILTER_QUERY.format(
            gene_array_filter=config.GENE_ARRAY)
    return db.iquery(
        query,
        schema=config.GENE_FILTER_SCHEMA,
        fetch=True,
        atts_only=True)


def get_gene_by_id(db, gene_id):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/gene/ENSG00000107404

    MongoDB:
      db.genes.find({'gene_id': 'ENSG00000107404'}, fields={'_id': False})

    SciDB:
      equi_join(
        equi_join(
            gene,
            filter(gene_index, gene_id = 'ENSG00000107404'),
            'left_names=gene_idx',
            'right_names=gene_idx',
            'algorithm=hash_replicate_right'),
        transcript_index,
        'left_names=transcript_idx',
        'right_names=transcript_idx',
        'algorithm=hash_replicate_right');
    """
    res = numpy2dict0(
        db.iquery(
            config.GENE_TRANSCRIPT_BY_ID_QUERY.format(gene_id=gene_id),
            schema=config.GENE_TRANSCRIPT_BY_ID_SCHEMA,
            fetch=True,
            atts_only=True))
    res['canonical_transcript'] = res['transcript_id']
    return res


def get_gene_by_id_new(db, gene_id):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/gene/ENSG00000107404

    MongoDB:
      db.genes.find({'gene_id': 'ENSG00000107404'}, fields={'_id': False})

    SciDB:
      equi_join(
        equi_join(
            gene,
            filter(gene_index, gene_id = 'ENSG00000107404'),
            'left_names=gene_idx',
            'right_names=gene_idx',
            'algorithm=hash_replicate_right'),
        transcript_index,
        'left_names=transcript_idx',
        'right_names=transcript_idx',
        'algorithm=hash_replicate_right');

    SciDBnew:
    genedf = bb.get_genes(gene_name='APOC3', exact_match=True)
    """
    genedf = db.get_genes(gene_name=gene_id, exact_match=True)
    return genedf


def exists_gene_id(db, gene_id):
    """
    Search bar

    MongoDB:
      db.genes.find({'gene_id': 'ENSG00000107404'}, fields={'_id': False})

    SciDB:
      aggregate(
        filter(gene_index, gene_id = 'ENSG00000198734'),
        count(*));

    SciDBnew: 
    ## needs to be updated
    res = bb.get_genes(gene_id='ENSG00000198734', exact_match=True)
    bool(res.empty)
    ### True or False 
    """
    res = db.get_genes(gene_id=gene_id, exact_match=True)
    return res

def get_genes_in_region(db, chrom, start, stop):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/region/16-50727514-50766988

    MongoDB:
      db.genes.find({'xstart': {'$lte': 1650766988},
                     'xstop' : {'$gte': 1650727514}}, fields={'_id': False})

    SciDB:
      equi_join(
        filter(gene,
               chrom = 16 and start <= 50766988 and stop >= 50727514),
        gene_index,
        'left_names=gene_idx',
        'right_names=gene_idx',
        'algorithm=hash_replicate_right');
    SciDBnew:
    bb.get_genes_in_region(namespace='RIVAS_HG38', chromosome=2, start=202201384, end=202201886)
    """
    res = db.get_genes_in_region(namespace=str(db.namespace), chromosome=int(chrom), start=int(start), end=int(stop))
    logging.info(res)
    return res

def get_gene_id_by_name(db, gene_name):
    """
    Search bar

    MondoDB:
      db.genes.find_one({'gene_name': 'F5'},   fields={'_id': False})
      db.genes.find_one({'other_names': 'F5'}, fields={'_id': False})

    SciDB:
      project(
        equi_join(
          filter(gene, gene_name = 'F5'),
          gene_index,
          'left_names=gene_idx',
          'right_names=gene_idx',
          'algorithm=hash_replicate_right'),
        gene_id);

    SciDBnew: 
     bb.get_genes(gene_name='RP11-150', exact_match=True)
    """
    # TODO - this needs to be updated 7.18.2019 MRivas
    res = db.get_genes(gene_name = gene_name, exact_match=True)
    if res.empty:
        return None
    return res['name'][0]


def exists_transcript_id(db, transcript_id):
    """
    Search bar

    MongoDB:
      db.transcripts.find({'transcript_id': 'ENST00000450546'},
                          fields={'_id': False})

    SciDB:
      aggregate(
        filter(transcript_index, transcript_id = 'ENST00000450546'),
        count(*));

    SciDBnew: 
    res = bb.get_transcript(transcript_id='RP11-150', exact_match=True)
    ## needs to be updated
    bool(res.empty)
    ### True or False 

    """
    res = db.get_transcripts(str(db.namespace), transcript_eid=str(transcript_id))
    if res.empty:
        return None
    else:
        return res

def get_transcript_gene(db, transcript_id):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/gene/ENSG00000107404

    MongoDB:
      db.transcripts.find({'transcript_id': 'ENST00000378891'},
                          fields={'_id': False})
      db.exons.find({'transcript_id': transcript_id,
                     'feature_type': { "$in": ['CDS', 'UTR', 'exon'] }},
                    fields={'_id': False})
      db.genes.find({'gene_id': 'ENSG00000107404'}, fields={'_id': False})
      db.transcripts.find({'gene_id': 'ENSG00000107404'},
                          fields={'_id': False})

    SciDB:
      See TRANSCRIPT_GENE_LOOKUP in config.py
    """
    return add_xpos(
        numpy2dict0(
            db.iquery(
                config.TRANSCRIPT_GENE_LOOKUP.format(
                    transcript_id=transcript_id),
                schema=config.TRANSCRIPT_GENE_SCHEMA,
                fetch=True,
                atts_only=True)))


# -- -
# -- - VARIANT - --
# -- -
def get_variants_chrom_pos_by_rsid_limit2(db, rsidentifier):
    """
    Search bar

    MongoDB:
      db.variants.find({'rsid': 'rs6025'}, fields={'_id': False}))

    SciDB:
          limit(
        project(
          filter(variant, rsid = 6025),
          rsid),
        2);
    SciDBnew:
     bb.get_variants(rsid=['rs34221567', 'rs148534902', 'rs200712517'])

    """
    if rsidentifier.startswith('rs'):
        res = db.get_variants(rsid = rsidentifier)
    else:
        return None
    if res.empty:
        return None
    return res


def get_variant_ann_by_chrom_pos(db, chrom, start, ref, alt):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/variant/1-39381448

    MongoDB:
      db.variants.find({'xpos': '1039381448'}, fields={'_id': False})

    SciDB:
      between(vairant, 1, 39381448,
                       1, 39381448);

    SciDBnew: 
    db.get_variants(chromosome=22, position=(32334104, 38877475))
    
    """
    # variants = format_variants(
    #     numpy2dict(
    #         db.iquery(
    #             config.VARIANT_LIMIT_QUERY.format(chrom=chrom, start=start),
    #             schema=config.VARIANT_LIMIT_SCHEMA,
    #             fetch=True)),
    #     add_ann=True)
    vl = pandas.DataFrame([[chrom, start, ref, alt]], columns = ['chrom','pos','ref','alt']) 
    if 'enrichlogor' in list(db.get_variant_fields()['name']):
        #variants = db.get_variants(chromosome=chrom, position=start, variant_fields = ['chrom','pos','ref','alt', 'variant_identity', 'major_consequence', 'category', 'gene_name', 'gene_symbol', 'HGVSp', 'HGVSc',  'consequence', 'all_filters',  'annotations', 'maf', 'ukbb_freq', 'ld', 'rsid', 'enrichlogor','enrichp'])
        variants = db.get_variants(variant_list = vl,  variant_fields = ['ref','alt', 'variant_identity', 'major_consequence', 'category', 'gene_name', 'gene_symbol', 'HGVSp', 'HGVSc',  'consequence', 'all_filters',  'annotations', 'maf', 'ukbb_freq', 'ld', 'rsid', 'enrichlogor','enrichp'])
    else:
#        variants = db.get_variants(chromosome=chrom, position=start, variant_fields = ['chrom','pos','ref','alt', 'variant_identity', 'major_consequence', 'category', 'gene_name', 'gene_symbol', 'HGVSp', 'HGVSc',  'consequence', 'all_filters',  'annotations', 'maf', 'ukbb_freq', 'ld', 'rsid'])
        variants = db.get_variants(variant_list = vl, variant_fields = ['ref','alt', 'variant_identity', 'major_consequence', 'category', 'gene_name', 'gene_symbol', 'HGVSp', 'HGVSc',  'consequence', 'all_filters',  'annotations', 'maf', 'ukbb_freq', 'ld', 'rsid'])
   # refl = [ref]
   # altl = [alt]
   # variants.query('ref == @refl and alt == @altl', inplace = True)
    variant = variants if len(variants) else None
    if variant is None:
        return variant
#    if variant['rsid'] == '.' or variant['rsid'] is None:
#        rsid = get_dbsnp(db, chrom, start)
#        if rsid:
#            variant['rsid'] = 'rs{}'.format(rsid)
    return variant


def get_variants_by_gene_idx(db, gene_idx, gene_id):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/gene/ENSG00000107404

    MongoDB:
      db.variants.find({'genes': 'ENSG00000107404'}, fields={'_id': False})

    SciDB:
      equi_join(
        variant,
        between(variant_gene,
                173, null,
                173, null),
        'left_names=chrom,pos',
        'right_names=chrom,pos');
    """
    return format_variants(
        numpy2dict(
            db.iquery(
                config.VARIANT_GENE_OR_TRANSCRIPT_LOOKUP.format(
                    variant_ref_array=config.VARIANT_GENE_ARRAY,
                    ref_idx=gene_idx),
                schema=config.VARIANT_GENE_OR_TRANSCRIPT_SCHEMA,
                fetch=True,
                atts_only=True)),
        gene_id=gene_id)


def get_variants_by_transcript_idx(db, transcript_idx, transcript_id):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/gene/ENSG00000107404

    MongoDB:
      db.variants.find({'transcripts': 'ENST00000378891'},
                       fields={'_id': False})

    SciDB:
      equi_join(
        variant,
        between(variant_transcript,
                3694, null,
                3694, null),
        'left_names=chrom,pos',
        'right_names=chrom,pos');
    """
    return format_variants(
        numpy2dict(
            db.iquery(
                config.VARIANT_GENE_OR_TRANSCRIPT_LOOKUP.format(
                    variant_ref_array=config.VARIANT_TRANSCRIPT_ARRAY,
                    ref_idx=transcript_idx),
                schema=config.VARIANT_GENE_OR_TRANSCRIPT_SCHEMA,
                fetch=True,
                atts_only=True)),
        transcript_id=transcript_id)


def get_variants_in_region(db, chrom, start, stop):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/region/16-50727514-50766988

    MongoDB:
      db.variants.find({'xpos': {'$lte': 1650766988},
                                 '$gte': 1650727514}}, fields={'_id': False})

    SciDB:
      between(variant, 16, 50727514,
                       16, 50766988);
    SciDBnew:
    bb.get_variants(chromosome=(21, 22), start=32300000, end=32400000)
    """
    # TODO add SEARCH_LIMIT
    if stop is None:
        stop = start
    res = db.get_variants(chromosome=int(chrom), start=int(start), end=int(stop), variant_fields = ('chrom','pos','ref','alt', 'variant_identity', 'major_consequence', 'category', 'gene_name', 'gene_symbol', 'HGVSp', 'HGVSc',  'consequence', 'all_filters',  'maf', 'ukbb_freq', 'ld', 'rsid'))
    return res

# -- -
# -- - COVERAGE - --
# -- -
def get_coverage_for_transcript(db, xstart, xstop=None):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/gene/ENSG00000107404

    MongoDB:
      db.base_coverage.find({'xpos': {'$gte': xstart, '$lte': xstop}},
                            fields={'_id': False})

    SciDB:
      between(coverage, 1, 1270607,
                        1, 1284543);

    """
    # TODO has_coverage filter
    return get_coverage_for_bases(db, xstart, xstop)


def get_coverage_for_bases(db, xstart, xstop=None):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/region/16-50727514-50766988

    MongoDB:
      db.base_coverage.find({'xpos': {'$gte': xstart, '$lte': xstop}},
                            fields={'_id': False})

    SciDB:
      between(coverage, 16, 50727514,
                        16, 50766988);
    """
    if xstop is None:
        xstop = xstart
    return numpy2dict(
        db.iquery(
            config.COVERAGE_LOOKUP_QUERY.format(
                chrom_start=int(xstart / config.XOFF),
                pos_start=int(xstart % config.XOFF),
                chrom_stop=int(xstop / config.XOFF),
                pos_stop=int(xstop % config.XOFF)),
            schema=config.COVERAGE_SCHEMA_OBJ,
            fetch=True))


# -- -
# -- - DBSNP - --
# -- -
def get_dbsnp(db, chrom, pos):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/variant/1-39381448

    MongoDB:
      db.dbsnp.find_one({'xpos': '1039381448'})

    SciDB:
      between(dbsnp_by_chrom_pos, 1, 39381448,
                                  1, 39381448),
    """
    res = db.iquery(
        config.DBSNP_LOOKUP_QUERY.format(chrom=chrom, pos=pos),
        schema=config.DBSNP_BY_CHROM_POS_SCHEMA_OBJ,
        fetch=True,
        atts_only=True)
    if not res:
        return None
    return res[0]['rsid']['val']


def get_variants_from_dbsnp(db, rsid):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/awesome?query=rs1171

    MongoDB:
      db.dbsnp.find_one({'rsid': rsid})
      db.variants.find({'xpos': {'$lte': position['xpos'],
                                 '$gte': position['xpos']}},
                       fields={'_id': False})

    SciDB:
      equi_join(
        filter(dbsnp_by_rsid, rsid = 1171),
        project(variant, rsid),
        'left_names=chrom,pos',
        'right_names=chrom,pos',
        'algorithm=hash_replicate_left');

    SciDBnew:
     bb.get_variants(rsid=['rs34221567', 'rs148534902', 'rs200712517'])
    """
    if not rsid.startswith('rs'):
        return None
    rsidarr = []
    rsidarr.append(rsid)
    res = db.get_variants(rsid=rsidarr)
    if res.empty:
        return None
    return '{}-{}-{}-{}'.format(res['chrom'][0], res['pos'][0], res['ref'][0], res['alt'][0])


# -- -
# -- - SEARCH BAR - --
# -- -
def get_awesomebar_suggestions(g, query):
    """This generates autocomplete suggestions when user query is the
    string that user types If it is the prefix for a gene, return list
    of gene names

    """
    regex = re.compile('^' + re.escape(query), re.IGNORECASE)
    results = [r for r in g.autocomplete_strings if regex.match(r)][:50]
    return results


def get_awesomebar_result(db, query):
    """Similar to the above, but this is after a user types enter We need to
    figure out what they meant - could be gene, variant, region, phenotype

    Return tuple of (datatype, identifier)
    Where datatype is one of 'gene', 'variant', or 'region'
    And identifier is one of:
    - ensembl ID for gene
    - variant ID string for variant (eg. 1-1000-A-T)
    - region ID string for region (eg. 1-1000-2000)

    Follow these steps:
    - if query is an ensembl ID, return it
    - if a gene symbol, return that gene's ensembl ID
    - if an RSID, return that variant's string
    - if a phenotype ID, return that phenotype

    Finally, note that we don't return the whole object here - only
    it's identifier.  This could be important for performance later

    """
    query = query.strip()
    (query_lower, query_upper) = (query.lower(), query.upper())
    query = str(query)
    query_lower = str(query_lower)
    query_upper = str(query_upper)
    if query_upper in UNSUPPORTED_QUERIES:
        return 'error', query

    # Variant
    variants = get_variants_chrom_pos_by_rsid_limit2(db, query)
    if variants is not None and not variants.empty:
        if len(variants) == 1:
            variant = variants            
            variant_id = '{}-{}-{}-{}'.format(variant['chrom'][0], variant['pos'][0], variant['ref'][0], variant['alt'][0])
            return 'variant', variant_id
        else:
            return 'dbsnp_variant_set', query_lower
    ## variant id
    variant_id = get_variants_from_dbsnp(db, query_lower)
    if variant_id:
        return 'variant', variant_id
    ## gene id - here now - needs to be updated
    gene_id = get_gene_id_by_name(db, query_upper)
    if gene_id:
        return 'gene', gene_id

    # Ensembl formatted queries - here now - needs to be updated
    if query_upper.startswith('ENS'):
        # Gene
        if query_upper.startswith('ENSG'):
            if exists_gene_id(db, query_upper):
                return 'gene', query_upper

        # Transcript - here now - needs to be updated
        if exists_transcript_id(db, query_upper) is not None:
            res = exists_transcript_id(db, query_upper) 
            return 'region', '{}-{}-{}'.format(res['chrom'][0], res['start'][0], res['end'][0])

    # GBE ID formatted queries - here now - needs to be updated
    if exists_icd(db, query):
        # HC382
        phe_title = get_phe_title(db, query) 
        return 'gbe', phe_title

    # From here on out, only region queries
    if query_upper.startswith('CHR'):
        query_upper = query_upper.lstrip('CHR')

    # Region
    m = REGION_RE1.match(query_upper)
    if m:
        if int(m.group(3)) < int(m.group(2)):
            return 'region', 'invalid'
        return 'region', '{}-{}-{}'.format(m.group(1), m.group(2), m.group(3))
    m = REGION_RE2.match(query_upper)
    if m:
        return 'region', '{}-{}-{}'.format(m.group(1), m.group(2), m.group(2))
    m = REGION_RE3.match(query_upper)
    if m:
        return 'region', '{}'.format(m.group(1))
    m = REGION_RE4.match(query_upper)
    if m:
        return 'variant', '{}-{}-{}-{}'.format(
            m.group(1), m.group(2), m.group(3), m.group(4))

    return 'not_found', query


# -- -
# -- - MODELS - --
# -- -
def get_gene_variant(db, gene_names=None, icds=None):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/runPolyCoding

    MongoDB:
      db.genes.find_one({'gene_name': gene_name}, fields={'_id': False})
      db.variants.find({'genes': gene_id}, fields={'_id': False})

    SciDB:
      See GENE_VARIANT_LOOKUP in config.py
    """
    if gene_names:
        gene_filter = 'filter({gene_array}, {cond})'.format(
            gene_array=config.GENE_ARRAY,
            cond=' or '.join(
                "gene_name = '{}'".format(g) for g in gene_names))
    else:
        gene_filter = config.GENE_ARRAY

    return db.iquery(
        config.GENE_VARIANT_LOOKUP.format(gene_filter=gene_filter, bim_filter="true", variant_filter="true"),
        #schema=config.GENE_VARIANT_SCHEMA,
        fetch=True,
        atts_only=True)


def get_variant_icd(db, gene_names=None, icds=None):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/runPolyCoding

    MongoDB:
      db.genes.find_one({'gene_name': gene_name}, fields={'_id': False})
      db.variants.find({'genes': gene_id}, fields={'_id': False})
      db.icd.find({'xpos': {'$gte': gene.xstart, '$lte': gene.xstop}},
                  fields={'_id': False})

    SciDB:
      See VARIANT_ICD_LOOKUP in config.py
    """
    if gene_names:
        gene_filter = 'filter({gene_array}, {cond})'.format(
            gene_array=config.GENE_ARRAY,
            cond=' or '.join(
                "gene_name = '{}'".format(g) for g in gene_names))
    else:
        gene_filter = config.GENE_ARRAY

    if icds:
        icd_filter = 'filter({icd_info_array}, {cond})'.format(
            icd_info_array=config.ICD_INFO_ARRAY,
            cond=' or '.join("icd = '{}'".format(i) for i in icds))
    else:
        icd_filter = config.ICD_INFO_ARRAY

    return db.iquery(
        config.VARIANT_ICD_LOOKUP.format(icd_filter=icd_filter,
                                         gene_filter=gene_filter),
        schema=config.VARIANT_ICD_SCHEMA,
        fetch=True,
        atts_only=True)


def print_all():
    import pprint

    pp = pprint.PrettyPrinter(indent=2)
    db = scidbpy.connect()

    pp.pprint(get_icd_name_map(db))

    # /coding/RH117 -> gbe.icd_page()
    pp.pprint(get_icd_variant_by_icd_id_pvalue(db, 'RH117'))

    # /variant/1-169519049 -> gbe.variant_icd_page()
    pp.pprint(get_variant_ann_by_chrom_pos(db, 1, 169519049))
    pp.pprint(get_icd_by_chrom_pos(db, 1, 169519049))

    # /gene/ENSG00000107404 -> gbe.gene_page()
    pp.pprint(get_gene_by_id(db, 'ENSG00000107404'))
    pp.pprint(get_variants_by_gene_idx(db, 173, 'ENSG00000107404'))
    pp.pprint(get_variants_by_transcript_idx(db, 3694, 'ENST00000378891'))
    pp.pprint(get_coverage_for_transcript(db, 1001270606, 1001284542))

    # /transcript/ENST00000289248 -> gbe.transcript_page()
    pp.pprint(get_transcript_gene(db, 'ENST00000289248'))
    pp.pprint(get_variants_by_transcript_idx(db, 304, 'ENST00000289248'))
    pp.pprint(get_coverage_for_transcript(db, 1039351918, 1039392559))

    # /region/1-28052491-28089634 -> gbe.region_page()
    pp.pprint(get_genes_in_region(db, 1, 28052491, 28089634))
    pp.pprint(get_variants_in_region(db, 1, 28052491, 28089634))
    pp.pprint(get_coverage_for_bases(db, 1028052491, 1028089634))

    # /region/16-50727514-50766988 -> gbe.region_page()
    pp.pprint(get_genes_in_region(db, 16, 50727514, 50766988))
    pp.pprint(get_variants_in_region(db, 16, 50727514, 50766988))
    pp.pprint(get_coverage_for_bases(db, 16050727514, 16050766988))

    # /intensity/ -> gbe.intensity_page()
    pp.pprint(get_icd_affyid(db, 723307))

    # /awesome -> gbe.awesome()
    pp.pprint(get_variants_chrom_pos_by_rsid_limit2(db, 'rs6025'))
    pp.pprint(get_variants_from_dbsnp(db, 'rs1171'))
    pp.pprint(get_gene_id_by_name(db, 'F5'))
    pp.pprint(exists_gene_id(db, 'ENSG00000107404'))
    pp.pprint(exists_transcript_id(db, 'ENST00000378891'))
    pp.pprint(exists_icd(db, 'RH117'))

    # /target/1 -> gbe.target_page()
    # pp.pprint(get_icd_variant_by_pvalue(db))

    # models
    pp.pprint(get_gene_variant(db))
    pp.pprint(get_variant_icd(db, icds=('RH117',)))


def time_all():
    import timer

    db = None

    #
    # init
    #
    tm = 0
    with timer.Timer() as t:
        db = scidbpy.connect()
    print('{:8.2f}ms\t{}'.format(t.msecs, 'connect'))
    tm += t.msecs

    with timer.Timer() as t:
        get_icd_name_map(db)
    print('{:8.2f}ms\t{}'.format(t.msecs, 'get_icd_name_map'))
    tm += t.msecs
    print('     -----\n{:8.2f}ms\t{}\n'.format(tm, 'init'))

    #
    # /coding/RH117 -> gbe.icd_page()
    #
    with timer.Timer() as t:
        get_icd_variant_by_icd_id_pvalue(db, 'RH117')
    print('{:8.2f}ms\t{}'.format(t.msecs, 'get_icd_variant_by_icd_id_pvalue'))
    print('     -----\n{:8.2f}ms\t{}\n'.format(t.msecs, '/coding/'))

    #
    # /variant/1-169519049 -> gbe.variant_icd_page()
    #
    tm = 0
    with timer.Timer() as t:
        get_variant_ann_by_chrom_pos(db, 1, 169519049)
    print('{:8.2f}ms\t{}'.format(t.msecs, 'get_variant_ann_by_chrom_pos'))
    tm += t.msecs

    with timer.Timer() as t:
        get_icd_by_chrom_pos(db, 1, 169519049)
    print('{:8.2f}ms\t{}'.format(t.msecs, 'get_icd_by_chrom_pos'))
    tm += t.msecs
    print('     -----\n{:8.2f}ms\t{}\n'.format(tm, '/variant/'))

    #
    # /gene/ENSG00000107404 -> gbe.gene_page()
    #
    tm = 0
    with timer.Timer() as t:
        get_gene_by_id(db, 'ENSG00000107404')
    print('{:8.2f}ms\t{}'.format(t.msecs, 'get_gene_by_id'))
    tm += t.msecs

    with timer.Timer() as t:
        get_variants_by_gene_idx(db, 173, 'ENSG00000107404')
    print('{:8.2f}ms\t{}'.format(t.msecs, 'get_variants_by_gene_idx'))
    tm += t.msecs

    with timer.Timer() as t:
        get_variants_by_transcript_idx(db, 3694, 'ENST00000378891')
    print('{:8.2f}ms\t{}'.format(t.msecs, 'get_variants_by_transcript_idx'))
    tm += t.msecs

    with timer.Timer() as t:
        get_coverage_for_transcript(db, 1001270607, 1001284543)
    print('{:8.2f}ms\t{}'.format(t.msecs, 'get_coverage_for_transcript'))
    tm += t.msecs
    print('     -----\n{:8.2f}ms\t{}\n'.format(tm, '/gene/'))

    #
    # /transcript/ENST00000289248 -> gbe.transcript_page()
    #
    tm = 0
    with timer.Timer() as t:
        get_transcript_gene(db, 'ENST00000289248')
    print('{:8.2f}ms\t{}'.format(t.msecs, 'get_transcript_gene'))
    tm += t.msecs

    with timer.Timer() as t:
        get_variants_by_transcript_idx(db, 304, 'ENST00000289248')
    print('{:8.2f}ms\t{}'.format(t.msecs, 'get_variants_by_transcript_idx'))
    tm += t.msecs

    with timer.Timer() as t:
        get_coverage_for_transcript(db, 1039351919, 1039392560)
    print('{:8.2f}ms\t{}'.format(t.msecs, 'get_coverage_for_transcript'))
    tm += t.msecs
    print('     -----\n{:8.2f}ms\t{}\n'.format(tm, '/transcript/'))

    #
    # /region/1-28052491-28089634 -> gbe.region_page()
    #
    tm = 0
    with timer.Timer() as t:
        get_genes_in_region(db, 1, 28052491, 28089634)
    print('{:8.2f}ms\t{}'.format(t.msecs, 'get_genes_in_region'))
    tm += t.msecs

    with timer.Timer() as t:
        get_variants_in_region(db, 1, 28052491, 28089634)
    print('{:8.2f}ms\t{}'.format(t.msecs, 'get_variants_in_region'))
    tm += t.msecs

    with timer.Timer() as t:
        get_coverage_for_bases(db, 1028052491, 1028089634)
    print('{:8.2f}ms\t{}'.format(t.msecs, 'get_coverage_for_bases'))
    tm += t.msecs
    print('     -----\n{:8.2f}ms\t{}\n'.format(tm, '/region/'))

    #
    # /intensity/ -> gbe.intensity_page()
    #
    with timer.Timer() as t:
        get_icd_affyid(db, 723307)
    print('{:8.2f}ms\t{}'.format(t.msecs, 'get_icd_affyid'))
    print('     -----\n{:8.2f}ms\t{}\n'.format(t.msecs, '/intensity/'))

    #
    # /awesome -> gbe.awesome()
    #
    tm = 0
    with timer.Timer() as t:
        get_variants_chrom_pos_by_rsid_limit2(db, 'rs6025')
    print('{:8.2f}ms\t{}'.format(t.msecs,
                                 'get_variants_chrom_pos_by_rsid_limit2'))
    tm += t.msecs

    with timer.Timer() as t:
        get_variants_from_dbsnp(db, 'rs771157073')
    print('{:8.2f}ms\t{}'.format(t.msecs, 'get_variants_from_dbsnp'))
    tm += t.msecs

    with timer.Timer() as t:
        get_gene_id_by_name(db, 'F5')
    print('{:8.2f}ms\t{}'.format(t.msecs, 'get_gene_id_by_name'))
    tm += t.msecs

    with timer.Timer() as t:
        exists_gene_id(db, 'ENSG00000107404')
    print('{:8.2f}ms\t{}'.format(t.msecs, 'exists_gene_id'))
    tm += t.msecs

    with timer.Timer() as t:
        exists_transcript_id(db, 'ENST00000378891')
    print('{:8.2f}ms\t{}'.format(t.msecs, 'exists_transcript_id'))
    tm += t.msecs

    with timer.Timer() as t:
        exists_icd(db, 'RH117')
    print('{:8.2f}ms\t{}'.format(t.msecs, 'exists_icd'))
    tm += t.msecs
    print('     -----\n{:8.2f}ms\t{}\n'.format(tm, '/awesome'))

    #
    # /target/1 -> gbe.target_page()
    #
    # with timer.Timer() as t:
    #     get_icd_variant_by_pvalue(db)
    # print('{:8.2f}ms\t{}'.format(t.msecs, 'get_icd_variant_by_pvalue'))
    # print('     -----\n{:8.2f}ms\t{}\n'.format(t.msecs, '/target/1'))

    #
    # /runPolyCoding -> gbe.runPolyCoding_page()
    #
    tm = 0
    with timer.Timer() as t:
        get_gene_variant(db)
    print('{:8.2f}ms\t{}'.format(t.msecs, 'get_gene_variant'))
    tm += t.msecs

    with timer.Timer() as t:
        get_variant_icd(db, icds=('RH117',))
    print('{:8.2f}ms\t{}'.format(t.msecs, 'get_variant_icd'))
    tm += t.msecs
    print('     -----\n{:8.2f}ms\t{}\n'.format(tm, '/runPolyCoding'))


if __name__ == '__main__':
    print_all()
