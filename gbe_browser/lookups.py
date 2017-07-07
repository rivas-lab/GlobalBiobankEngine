import itertools
import scidbpy


import config
import utils


def numpy2dict(ar):
    """Convert SciDB NumPy array result to Python dictionary and populate
    nullable attributes with values (discards null codes).

    """
    return [
        dict(
            (de[0],
             el[de[0]]['val'] if isinstance(de[1], list) else el[de[0]])
            for de in ar.dtype.descr)
        for el in ar]


def format_variants(variants, add_ann=False):
    for variant in variants:
        variant['rsid'] = 'rs{}'.format(variant['rsid'])
        variant['variant_id'] = '{}-{}-{}-{}'.format(
            variant['chrom'], variant['pos'], variant['ref'], variant['alt'])

        anns = [dict(zip(config.VARIANT_CSQ, csq.split('|')))
                for csq in variant['csq'].split(',')]
        vep_annotations = [ann for ann in anns
                           if ('Feature' in ann and
                               ann['Feature'].startswith('ENST'))]
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


# -- -
# -- - ICD - --
# -- -
def get_icd_significant(db, icd_id, cutoff=0.01):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/coding/RH117

    MongoDB:
      db.icd.find({'icd': 'RH117', 'stats.pvalue': {'$lt': 0.01}},
                  fields={'_id': false})

    SciDB:
      filter(
        cross_join(icd,
                   filter(icd_index, icd = 'RH117'),
                   icd.icd_idx,
                   icd_index.icd_idx),
        pvalue < 0.01);
    """
    return numpy2dict(
        db.iquery(
            config.ICD_PVALUE_LOOKUP_QUERY.format(
                icd=icd_id, pvalue=cutoff),
            schema=config.ICD_X_INFO_SCHEMA,
            fetch=True))


def get_icd_info(db, icd_id):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/coding/RH117

    MongoDB:
      db.icd_info.find({'icd': 'RH117'}, fields={'_id': False})

    SciDB:
      cross_join(icd_info,
                 filter(icd_index, icd = 'RH117'),
                 icd_info.icd_idx,
                 icd_index.icd_idx);
    """
    return numpy2dict(
        db.iquery(
            config.ICD_LOOKUP_QUERY.format(icd=icd_id),
            schema=config.ICD_X_INFO_SCHEMA,
            fetch=True))


def get_variant_icd(db, xpos):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/variant/1-39381448

    MongoDB:
      db.icd.find({'xpos': '1039381448'}, fields={'_id': False})

    SciDB:
      between(icd, null, 1, 39381448, null, null, 1, 39381448, null);
    """
    return numpy2dict(
        db.iquery(
            config.ICD_CHROM_POS_LOOKUP_QUERY.format(
                chrom=int(xpos / 1e9), pos=int(xpos % 1e9)),
            schema=config.ICD_X_INFO_SCHEMA,
            fetch=True))


def get_icd_significant_variant(db, icd_id, cutoff=0.001):
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
      cross_join(
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
        variant.chrom,
        icd_join.chrom,
        variant.pos,
        icd_join.pos);
    """
    pdecimal = config.ICD_PVALUE_MAP.get(cutoff, 0)
    return format_variants(numpy2dict(db.iquery(
        config.ICD_VARIANT_LOOKUP_QUERY.format(
            icd=icd_id, pdecimal=pdecimal),
        schema=config.VARIANT_X_ICD_X_INFO_SCHEMA,
        fetch=True)))


# -- -
# -- - VARIANT - --
# -- -
def get_variants_by_id(db, variant_ids):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/coding/RH117

    MongoDB:
      db.variants.find({'xpos': '1039381448'}, fields={'_id': False})

    SciDB:
      filter(variant, chrom = 1 and pos = 39381448); -- or chrom = ...
    """
    chrom_pos_cond = ' or '.join(
        'chrom = {chrom} and pos = {pos}'.format(
            chrom=int(xpos / 1e9), pos=int(xpos % 1e9))
        for xpos in variant_ids)
    variants = numpy2dict(
        db.iquery(
            config.VARIANT_MULTI_LOOKUP_QUERY.format(
                chrom_pos_cond=chrom_pos_cond),
            schema=config.VARIANT_LOOKUP_SCHEMA,
            fetch=True))
    variants = format_variants(variants)
    return variants


def get_variant_chrom_pos(db, chrom, pos):
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
    variants = numpy2dict(
        db.iquery(
            config.VARIANT_LOOKUP_QUERY.format(chrom=chrom, pos=pos),
            schema=config.VARIANT_LOOKUP_SCHEMA,
            fetch=True))
    variants = format_variants(variants, add_ann=True)
    variant = variants[0] if len(variants) else None
    if variant is None or 'rsid' not in variant:
        return variant
    if variant['rsid'] == '.' or variant['rsid'] is None:
        raise NotImplementedError()  # TODO
        # rsid = db.dbsnp.find_one({'xpos': xpos})
        # if rsid:
        #     variant['rsid'] = 'rs%s' % rsid['rsid']
    return variant


def get_variant(db, xpos):
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
    return get_variant_chrom_pos(db, int(xpos / 1e9), int(xpos % 1e9))


# -- -
# -- - GENE - --
# -- -
def get_gene(db, gene_id):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/gene/ENSG00000107404

    MongoDB:
      db.genes.find({'gene_id': 'ENSG00000107404'}, fields={'_id': False})

    SciDB:
      cross_join(gene,
                 filter(gene_index, gene_id = 'ENSG00000107404'),
                 gene.gene_idx,
                 gene_index.gene_idx);

    """
    return numpy2dict(
        db.iquery(
            config.GENE_LOOKUP_QUERY.format(gene_id=gene_id),
            schema=config.GENE_LOOKUP_SCHEMA,
            fetch=True))


def get_variants_in_gene(db, gene_id):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/gene/ENSG00000107404

    MongoDB:
      db.variants.find({'genes': 'ENSG00000107404'}, fields={'_id': False})

    SciDB:
      cross_join(variant,
                 cross_join(variant_gene,
                            filter(gene_index, gene_id = 'ENSG00000107404'),
                            variant_gene.gene_idx,
                            gene_index.gene_idx) as variant_gene_index,
                 variant.chrom,
                 variant_gene_index.chrom,
                 variant.pos,
                 variant_gene_index.pos);
    """
    return format_variants(
        numpy2dict(
            db.iquery(
                config.VARIANT_GENE_LOOKUP.format(gene_id=gene_id),
                schema=config.VARIANT_X_GENE_INDEX_SCHEMA,
                fetch=True)))


if __name__ == '__main__':
    db = scidbpy.connect()

    import pprint
    pp = pprint.PrettyPrinter(indent=2)
    pp.pprint(get_icd_significant(db, 'RH117'))
    pp.pprint(get_icd_info(db, 'RH117'))
    pp.pprint(get_variant_icd(db, 1039381448))
    pp.pprint(get_icd_significant_variant(db, 'RH117'))
    pp.pprint(get_variants_by_id(db, (1039381448,)))
    pp.pprint(get_variant(db, 1039381448))
    pp.pprint(get_gene(db, 'ENSG00000107404'))
    pp.pprint(get_variants_in_gene(db, 'ENSG00000107404'))
