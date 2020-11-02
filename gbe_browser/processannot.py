from __future__ import division

import lookups
import re
import scidbbiobank 
import config
import utils
import logging
import numpy 
import sys
import os
import scipy.stats

ns = sys.argv[1]

def get_db(name_space):
    DB = scidbbiobank.connect(scidb_url=os.getenv('SCIDB_URL',None), scidb_auth=('scidbadmin', 'Paradigm4'), namespace=name_space)
    DB.set_limit(15000)
    return DB

db = get_db(ns)


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


#def get_icd_variant_by_icd_id_pvalue(db, icd_id, field_identifier, pvalue=0.001):
assocset = str(db.list_association_sets()['name'][0])
chroms = range(1,23)
chroms.append('X')
chroms.append('Y')
for chrom in chroms:
    fnout = ns + '.' + str(chrom) + '.txt'
    """
    e.g.,OB
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
    df = db.get_variants(association_set=assocset, chromosome=chrom, variant_fields = ('chrom','pos','ref','alt', 'gnomad_filter',  'consequence', 'all_filters',  'annotations', 'maf', 'ld', 'rsid') )
    df['variant_id'] = df.apply(lambda row: RSID_FORMAT.format(
            chrom=row['chrom'],
            pos=row['pos'],
            ref=row['ref'],
            alt=row['alt']), axis = 1)
    df.loc[~df.consequence.isin(utils.csq_order),'consequence'] = 'intergenic'
    df['major_consequence'] = df['consequence']
    df['category'] = df.apply(lambda row: utils.add_category_to_variant(row['consequence']), axis = 1)
   # icdres['filter'] = icdres.apply(lambda row: 'PASS' if row['all_filters'] == 0 and row['se'] <= .5 else 'SE' if row['se'] > .5 else 'FAIL', axis = 1)
   # icdres['pvalue'] = icdres.apply(lambda row: row['pvalue'] if row['pvalue'] != 0 else 1e-300, axis = 1)
   # icdres['Name'] =  icdres.apply(lambda row: df[df['title'] == row['title']]['description'].squeeze(), axis = 1)
   # icdres['shortname'] = icdres.apply(lambda row: str(df[df['title'] == row['title']]['notes'].squeeze()).split(';')[1].split('=')[1], axis = 1)
   # icdres['Name'] = icdres['shortname']
   # icdres['Case'] = icdres.apply(lambda row: str(df[df['title'] == row['title']]['notes'].squeeze()).split(';')[0].split('=')[1], axis = 1)
    df['gene_name'], df['gene_symbol'], df['HGVSp'], df['HGVSc'] = zip(*df['annotations'].map(utils.return_gene_vep)) 
   # icdres = icdres.drop(columns=['annotations'])
    #icdres.apply(lambda row: utils.return_gene_vep(row['annotations']), axis = 1)
    #icdres['log10pvalue'] = icdres.apply(lambda row: -numpy.log10(row['pvalue']), axis = 1)
    #icdres['icd'] = icdres['title']
    #if icdres.loc[icdres['odds_ratio'].isna()].shape[0] < icdres.loc[icdres['beta'].isna()].shape[0]:
    #    icdres['or_val'] = icdres['odds_ratio']
    #    icdres['lor_val'] = icdres.apply(lambda row: numpy.log(row['odds_ratio']), axis = 1)
    #else:
    #    icdres['or_val'] = icdres['beta']
    #    icdres['lor_val'] = icdres['beta']
    df['ukbb_freq'] = df.apply(lambda row: min(float(row['maf']), 1 - float(row['maf'])) if row['maf'] is not None else None, axis = 1)
    df['gnomad_af'] = df.apply(lambda row: min(float(row['gnomad_filter']),1-float(row['gnomad_filter'])) if row['gnomad_filter'] is not None else None, axis = 1)
#    df['enrichlogor'] = df.apply(lambda row: numpy.log((float(row['ukbb_freq'])*(1-row['gnomad_af']))/((.00000000001 + row['gnomad_af'])*(0.00000000001 + 1-float(row['ukbb_freq'])))) if row['gnomad_af'] is not None else 10, axis = 1)
#    df['enrichp'] = df.apply(lambda row: scipy.stats.norm.sf(row['enrichlogor']/numpy.sqrt(1/(float(row['ukbb_freq'])*200000*2 + .5) + 1/((1-float(row['ukbb_freq']))*200000*2 + .5) + 1/(float(row['gnomad_af'])*30000*2 + .5) + 1/((1 - float(row['gnomad_af']))*30000*2 + .5))) if row['gnomad_filter'] is not None else 0, axis = 1)
  #  dfn = df[['chrom','pos', 'ref','alt','variant_id','major_consequence','category','gene_name','gene_symbol','HGVSp','HGVSc','ukbb_freq','gnomad_af','enrichlogor', 'enrichp']]
    dfn = df[['chrom','pos', 'ref','alt','variant_id','major_consequence','category','gene_name','gene_symbol','HGVSp','HGVSc','ukbb_freq','gnomad_af']]
    dfn.to_csv(fnout, index=None, sep='\t', mode='w')
