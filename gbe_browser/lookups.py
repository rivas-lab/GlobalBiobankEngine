import itertools
import re
import scidbpy

import config
import utils


UNSUPPORTED_QUERIES = set((
    'CMD1G',
    'CMH9',
    'CMPD4',
    'ENSG00000155657',
    'ENST00000342175',
    'ENST00000342992',
    'ENST00000359218',
    'ENST00000460472',
    'ENST00000589042',
    'ENST00000591111',
    'FLJ32040',
    'LGMD2J',
    'MYLK5',
    'TMD',
    'TTN',
))

RSID_FORMAT = '{chrom}-{pos}-{ref}-{alt}'

# 1:1-1000
REGION_RE1 = re.compile(r'^(\d+|X|Y|M|MT)\s*:\s*(\d+)-(\d+)$')
REGION_RE2 = re.compile(r'^(\d+|X|Y|M|MT)\s*:\s*(\d+)$')
REGION_RE3 = re.compile(r'^(\d+|X|Y|M|MT)$')
REGION_RE4 = re.compile(r'^(\d+|X|Y|M|MT)\s*[-:]\s*(\d+)-([ATCG]+)-([ATCG]+)$')


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


def format_gene(gene):
    gene['xstart'] = gene['chrom'] * config.XOFF + gene['start']
    gene['xstop'] = gene['chrom'] * config.XOFF + gene['stop']
    return gene


def format_genes(genes):
    for gene in genes:
        format_gene(gene)
    return genes


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
    """
    return exists(db, config.ICD_INFO_ARRAY, 'icd', "'{}'".format(icd))


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
    """
    if not stop:
        stop = start
    if not icd:
        icd_info_filter = config.ICD_INFO_ARRAY
    else:
        icd_info_filter = 'filter({icd_info_array}, {cond})'.format(
            icd_info_array=config.ICD_INFO_ARRAY,
            cond=' or '.join("icd = '{}'".format(i) for i in icd))
    return numpy2dict(
        db.iquery(
            config.ICD_CHROM_POS_LOOKUP_QUERY.format(
                icd_info_filter=icd_info_filter,
                chrom=chrom,
                start=start,
                stop=stop),
            schema=config.ICD_CHROM_POS_LOOKUP_SCHEMA,
            fetch=True,
            atts_only=True))


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


def get_icd_variant_by_icd_id_pvalue(db, icd_id, pvalue=0.001):
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
                config.ICD_VARIANT_LOOKUP_QUERY.format(
                    icd=icd_id, pdecimal=pdecimal),
                schema=config.VARIANT_X_ICD_X_INFO_SCHEMA,
                fetch=True,
                atts_only=True)))


def get_icd_variant_by_pvalue(db, pvalue=0.001):
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
            filter(icd_info, Cases >= 100),
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
                config.ICD_VARIANT_SCAN_QUERY.format(pdecimal=pdecimal),
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


def exists_gene_id(db, gene_id):
    """
    Search bar

    MongoDB:
      db.genes.find({'gene_id': 'ENSG00000107404'}, fields={'_id': False})

    SciDB:
      aggregate(
        filter(gene_index, gene_id = 'ENSG00000198734'),
        count(*));
    """
    return exists(db,
                  config.GENE_INDEX_ARRAY,
                  'gene_id',
                  "'{}'".format(gene_id))


def get_gene_by_idx(db, gene_idx):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/transcript/ENST00000458200

    MongoDB:
      db.genes.find({'gene_id': 'ENSG00000107404'}, fields={'_id': False})

    SciDB:
      equi_join(
        between(gene, 173, 173),
        transcript_index,
        'left_names=transcript_idx',
        'right_names=transcript_idx',
        'algorithm=hash_replicate_right');
    """
    res = numpy2dict0(
        db.iquery(
            config.GENE_TRANSCRIPT_BY_IDX_QUERY.format(gene_idx=gene_idx),
            schema=config.GENE_TRANSCRIPT_BY_IDX_SCHEMA,
            fetch=True,
            atts_only=True))
    res['canonical_transcript'] = res['transcript_id']
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
    """
    return numpy2dict(
        db.iquery(
            config.GENE_REGION_QUERY.format(
                chrom=chrom, start=start, stop=stop),
            schema=config.GENE_REGION_SCHEMA,
            fetch=True,
            atts_only=True))


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
    """
    # TODO other_names
    res = db.iquery(
        config.GENE_ID_BY_NAME_QUERY.format(gene_name=gene_name),
        schema=config.GENE_ID_BY_NAME_SCHEMA,
        fetch=True,
        atts_only=True)
    if not res:
        return None
    return res[0]['gene_id']['val']


def get_transcript_by_idx(db, transcript_idx):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/gene/ENSG00000107404

    MongoDB:
      db.transcripts.find({'transcript_id': 'ENST00000378891'},
                          fields={'_id': False})

    SciDB:
      between(transcript, null, 3694,
                          null, 3694);
    """
    return format_gene(
        numpy2dict0(
            db.iquery(
                config.TRANSCRIPT_IDX_LOOKUP.format(
                    gene_idx='null',
                    transcript_idx=transcript_idx),
                schema=config.TRANSCRIPT_SCHEMA_OBJ,
                fetch=True)))


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
    """
    return exists(db,
                  config.TRANSCRIPT_INDEX_ARRAY,
                  'transcript_id',
                  "'{}'".format(transcript_id))


def get_transcript_gene(db, transcript_id):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/gene/ENSG00000107404

    MongoDB:
      db.transcripts.find({'transcript_id': 'ENST00000378891'},
                          fields={'_id': False})

    SciDB:
      cross_join(
        cross_join(
          transcript,
          filter(transcript_index, transcript_id = 'ENST00000378891'),
          transcript.transcript_idx,
          transcript_index.transcript_idx),
        gene_index,
        transcript.gene_idx,
        gene_index.gene_idx);
    """
    res = format_gene(
        numpy2dict0(
            db.iquery(
                config.TRANSCRIPT_GENE_LOOKUP.format(
                    transcript_id=transcript_id),
                schema=config.TRANSCRIPT_GENE_SCHEMA,
                fetch=True)))
    return res


def get_transcripts_by_gene_idx(db, gene_idx):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/gene/ENSG00000107404

    MongoDB:
      db.transcripts.find({'gene_id': 'ENSG00000107404'},
                          fields={'_id': False})

    SciDB:
      between(transcript, 173, null,
                          173, null);
    """
    return numpy2dict(
        db.iquery(
            config.TRANSCRIPT_IDX_LOOKUP.format(
                gene_idx=gene_idx,
                transcript_idx='null'),
            schema=config.TRANSCRIPT_SCHEMA_OBJ,
            fetch=True))


def get_transcripts_id_by_gene_idx(db, gene_idx):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/gene/ENSG00000107404

    MongoDB:
      db.transcripts.find({'gene_id': 'ENSG00000107404'},
                          fields={'_id': False})

    SciDB:
      cross_join(
        between(transcript, 173, null,
                            173, null),
        transcript_index,
        transcript.transcript_idx,
        transcript_index.transcript_idx);
    """
    return numpy2dict(
        db.iquery(
            config.TRANSCRIPT_ID_LOOKUP.format(gene_idx=gene_idx),
            schema=config.TRANSCRIPT_ID_SCHEMA,
            fetch=True))


def get_exons(db, transcript_idx):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/gene/ENSG00000107404

    MongoDB:
      db.exons.find({'transcript_id': transcript_id,
                     'feature_type': { "$in": ['CDS', 'UTR', 'exon'] }},
                    fields={'_id': False})

    SciDB:
      between(exon, null, 3694, null, null, null, null,
                    null, 3694, null, null, null, null);
    """
    return numpy2dict(
        db.iquery(
            config.TRANSCRIPT_OR_EXON_BETWEEN_QUERY.format(
                array_name=config.EXON_ARRAY,
                gene_idx='null',
                transcript_idx=transcript_idx),
            schema=config.EXON_SCHEMA_OBJ,
            fetch=True))


# -- -
# -- - VARIANT - --
# -- -
def get_variants_chrom_pos_by_rsid_limit2(db, rsid):
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
    """
    if not rsid.startswith('rs'):
        return None
    rsid_int = None
    try:
        rsid_int = int(rsid.lstrip('rs'))
    except Exception:
        return None

    res = db.iquery(
        config.VARIANT_CHROM_POS_BY_RSID_QUERY.format(rsid=rsid_int),
        schema=config.VARIANT_CHROM_POS_BY_RSID_SCHEMA,
        fetch=True)
    if not res:
        return None

    return res


def get_variant_ann_by_chrom_pos(db, chrom, start):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/variant/1-39381448

    MongoDB:
      db.variants.find({'xpos': '1039381448'}, fields={'_id': False})

    SciDB:
      between(vairant, 1, 39381448,
                       1, 39381448);
    """
    variants = format_variants(
        numpy2dict(
            db.iquery(
                config.VARIANT_LIMIT_QUERY.format(chrom=chrom, start=start),
                schema=config.VARIANT_LIMIT_SCHEMA,
                fetch=True)),
        add_ann=True)
    variant = variants[0] if len(variants) else None
    if variant is None or 'rsid' not in variant:
        return variant
    if variant['rsid'] == '.' or variant['rsid'] is None:
        rsid = get_dbsnp(db, chrom, start)
        if rsid:
            variant['rsid'] = 'rs{}'.format(rsid)
    return variant


def get_variants_by_gene_idx(db, gene_idx, gene_id):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/gene/ENSG00000107404

    MongoDB:
      db.variants.find({'genes': 'ENSG00000107404'}, fields={'_id': False})

    SciDB:
      cross_join(variant,
                 between(variant_gene, null, null, 173,
                                       null, null, 173),
                 variant.chrom,
                 variant_gene.chrom,
                 variant.pos,
                 variant_gene.pos);
    """
    return format_variants(
        numpy2dict(
            db.iquery(
                config.VARIANT_GENE_LOOKUP.format(gene_idx=gene_idx),
                schema=config.VARIANT_GENE_SCHEMA,
                fetch=True)),
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
      cross_join(
        variant,
        between(variant_transcript, null, null, 3694,
                                    null, null, 3694),
        variant.chrom,
        variant_transcript.chrom,
        variant.pos,
        variant_transcript.pos);
    """
    return format_variants(
        numpy2dict(
            db.iquery(
                config.VARIANT_TRANSCRIPT_IDX_LOOKUP.format(
                    transcript_idx=transcript_idx),
                schema=config.VARIANT_X_TRANSCRIPT_SCHEMA,
                fetch=True)),
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
    """
    # TODO add SEARCH_LIMIT
    if stop is None:
        stop = start
    return format_variants(
        numpy2dict(
            db.iquery(
                config.VARIANT_LOOKUP_QUERY.format(
                    chrom=chrom, start=start, stop=stop),
                schema=config.VARIANT_LOOKUP_SCHEMA,
                fetch=True)),
        add_ann=True)


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
      https://biobankengine.stanford.edu/awesome?query=rs771157073

    MongoDB:
      db.dbsnp.find_one({'rsid': rsid})
      db.variants.find({'xpos': {'$lte': position['xpos'],
                                 '$gte': position['xpos']}},
                       fields={'_id': False})

    SciDB:
      equi_join(
        filter(dbsnp_by_rsid, rsid = 771157073),
        project(variant, rsid),
        'left_names=chrom,pos',
        'right_names=chrom,pos',
        'algorithm=hash_replicate_left');
    """
    if not rsid.startswith('rs'):
        return None
    rsid_int = None
    try:
        rsid_int = int(rsid.lstrip('rs'))
    except Exception:
        return None

    res = db.iquery(
        config.DBSNP_VARIANT_LOOKUP_QUERY.format(rsid=rsid_int),
        schema=config.DBSNP_VARIANT_LOOKUP_SCHEMA,
        fetch=True,
        atts_only=True)
    if not res:
        return None
    res0 = res[0]
    return '{}-{}'.format(res0['chrom']['val'], res0['pos']['val'])


# -- -
# -- - SEARCH BAR - --
# -- -
def get_awesomebar_suggestions(g, query):
    """This generates autocomplete suggestions when user query is the
    string that user types If it is the prefix for a gene, return list
    of gene names

    """
    regex = re.compile('^' + re.escape(query), re.IGNORECASE)
    results = [r for r in g.autocomplete_strings if regex.match(r)][:20]
    return results


def get_awesomebar_result(db, query):
    """Similar to the above, but this is after a user types enter We need to
    figure out what they meant - could be gene, variant, region, icd10

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


    Finally, note that we don't return the whole object here - only
    it's identifier.  This could be important for performance later

    """
    query = query.strip()
    (query_lower, query_upper) = (query.lower(), query.upper())
    print 'Query: %s' % query

    if query_upper in UNSUPPORTED_QUERIES:
        return 'error', query

    # Variant
    variants = get_variants_chrom_pos_by_rsid_limit2(db, query_lower)
    if variants:
        if len(variants) == 1:
            variant = variants[0]
            variant_id = '{}-{}'.format(variant['chrom'], variant['pos'])
            return 'variant', variant_id
        else:
            return 'dbsnp_variant_set', query_lower

    variant_id = get_variants_from_dbsnp(db, query.lower())
    if variant_id:
        return 'variant', variant_id

    gene_id = get_gene_id_by_name(db, query_upper)
    if gene_id:
        return 'gene', gene_id

    # Ensembl formatted queries
    if query_upper.startswith('ENS'):
        # Gene
        if exists_gene_id(db, query_upper):
            return 'gene', query_upper

        # Transcript
        if exists_transcript_id(db, query_upper):
            return 'transcript', query_upper

    # ICD10 formatted queries
    if query_upper.startswith('ICD'):
        # ICD10
        if exists_icd(db, query_upper):
            return 'icd10', query_upper

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
        config.GENE_VARIANT_LOOKUP.format(gene_filter=gene_filter),
        schema=config.GENE_VARIANT_SCHEMA,
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
    pp.pprint(get_transcripts_by_gene_idx(db, 173))
    pp.pprint(get_transcript_by_idx(db, 3694))
    pp.pprint(get_exons(db, 3694))
    pp.pprint(get_variants_by_transcript_idx(db, 3694, 'ENST00000378891'))
    pp.pprint(get_coverage_for_transcript(db, 1001270607, 1001284543))

    # /transcript/ENST00000289248 -> gbe.transcript_page()
    pp.pprint(get_transcript_gene(db, 'ENST00000289248'))
    pp.pprint(get_exons(db, 304))
    pp.pprint(get_gene_by_idx(db, 842))
    pp.pprint(get_transcripts_id_by_gene_idx(db, 842))
    pp.pprint(get_variants_by_transcript_idx(db, 304, 'ENST00000289248'))
    pp.pprint(get_coverage_for_transcript(db, 1039351919, 1039392560))

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
    pp.pprint(get_variants_from_dbsnp(db, 'rs771157073'))
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
        get_transcripts_by_gene_idx(db, 173)
    print('{:8.2f}ms\t{}'.format(t.msecs, 'get_transcripts_by_gene_idx'))
    tm += t.msecs

    with timer.Timer() as t:
        get_transcript_by_idx(db, 3694)
    print('{:8.2f}ms\t{}'.format(t.msecs, 'get_transcript_by_idx'))
    tm += t.msecs

    with timer.Timer() as t:
        get_exons(db, 3694)
    print('{:8.2f}ms\t{}'.format(t.msecs, 'get_exons'))
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
        get_exons(db, 304)
    print('{:8.2f}ms\t{}'.format(t.msecs, 'get_exons'))
    tm += t.msecs

    with timer.Timer() as t:
        get_gene_by_idx(db, 842)
    print('{:8.2f}ms\t{}'.format(t.msecs, 'get_gene_by_idx'))
    tm += t.msecs

    with timer.Timer() as t:
        get_transcripts_id_by_gene_idx(db, 842)
    print('{:8.2f}ms\t{}'.format(t.msecs, 'get_transcripts_id_by_gene_idx'))
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
    with timer.Timer() as t:
        get_icd_variant_by_pvalue(db)
    print('{:8.2f}ms\t{}'.format(t.msecs, 'get_icd_variant_by_pvalue'))
    print('     -----\n{:8.2f}ms\t{}\n'.format(t.msecs, '/target/1'))

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
