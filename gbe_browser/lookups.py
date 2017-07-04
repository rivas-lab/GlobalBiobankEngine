import exceptions
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
            schema=config.ICD_LOOKUP_SCHEMA_INST,
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
            config.ICD_INFO_LOOKUP_QUERY.format(icd=icd_id),
            schema=config.ICD_INFO_LOOKUP_SCHEMA_INST,
            fetch=True))


def get_variant_icd(db, xpos):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/variant/1-39381448

    MongoDB:
      db.icd.find({'xpos': '1039381448'}, fields={'_id': False})

    SciDB:
      filter(icd, xpos = 1039381448);
    """
    return numpy2dict(
        db.iquery(
            config.ICD_XPOS_LOOKUP_QUERY.format(xpos=xpos),
            schema=config.ICD_LOOKUP_SCHEMA_INST,
            fetch=True))


# -- -
# -- - VARIANT - --
# -- -
def format_variants(variants):
    for variant in variants:
        variant['rsid'] = 'rs{}'.format(variant['rsid'])
        variant['variant_id'] = '{}-{}-{}-{}'.format(
            variant['chrom'], variant['pos'], variant['ref'], variant['alt'])

        anns = [dict(zip(config.VARIANT_CSQ, csq.split('|')))
                for csq in variant['csq'].split(',')]
        vep_annotations = [ann for ann in anns
                           if ('Feature' in ann and
                               ann['Feature'].startswith('ENST'))]
        variant['vep_annotations'] = vep_annotations
        variant['genes'] = list(set(ann['Gene'] for ann in vep_annotations))
        variant['transcripts'] = list(set(
            ann['Feature'] for ann in vep_annotations))
    return variants


def get_variants_by_id(db, variant_ids):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/coding/RH117

    MongoDB:
      db.variants.find({'xpos': '1039381448'}, fields={'_id': False})

    SciDB:
      filter(variant, xpos = 1039381448);
    """
    xpos_cond = ' or '.join('xpos = {}'.format(xpos)
                            for xpos in variant_ids)
    variants = numpy2dict(
        db.iquery(
            config.VARIANT_LOOKUP_QUERY.format(xpos_cond=xpos_cond),
            schema=config.VARIANT_LOOKUP_SCHEMA_INST,
            fetch=True))
    variants = format_variants(variants)
    for variant in variants:
        utils.add_consequence_to_variant(variant)
    return variants


def get_variant(db, xpos):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/variant/1-39381448

    MongoDB:
      db.variants.find({'xpos': '1039381448'}, fields={'_id': False})

    SciDB:
      cross_join(filter(icd, xpos = 1039381448),
                 icd_index,
                 icd_info.icd_idx,
                 icd_index.icd_idx);
    """
    xpos_cond = 'xpos = {}'.format(xpos)
    variants = numpy2dict(
        db.iquery(
            config.VARIANT_LOOKUP_QUERY.format(xpos_cond=xpos_cond),
            schema=config.VARIANT_LOOKUP_SCHEMA_INST,
            fetch=True))
    variants = format_variants(variants)
    variant = variants[0] if len(variants) else None
    if variant is None or 'rsid' not in variant:
        return variant
    if variant['rsid'] == '.' or variant['rsid'] is None:
        raise NotImplementedError()  # TODO
        # rsid = db.dbsnp.find_one({'xpos': xpos})
        # if rsid:
        #     variant['rsid'] = 'rs%s' % rsid['rsid']
    return variant


# -- -
# -- - GENE - --
# -- -
def get_gene(db, gene_id):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/gene/ENSG00000101255

    MongoDB:
      db.genes.find({'gene_id': 'ENSG00000101255'}, fields={'_id': False})

    SciDB:
      cross_join(gene,
                 filter(gene_index, gene_id = 'ENSG00000101255'),
                 gene.gene_idx,
                 gene_index.gene_idx);

    """
    return numpy2dict(
        db.iquery(
            config.GENE_LOOKUP_QUERY.format(gene_id=gene_id),
            schema=config.GENE_LOOKUP_SCHEMA_INST,
            fetch=True))


def get_variants_in_gene(db, gene_id):
    """
    e.g.,
    UI:
      https://biobankengine.stanford.edu/gene/ENSG00000101255

    MongoDB:
      db.variants.find({'genes': 'ENSG00000101255'}, fields={'_id': False})

    SciDB:
      TODO
    """
    # TODO
    # variants = []
    # for variant in db.variants.find({'genes': gene_id},
    #                                 fields={'_id': False}):
    #     variant['vep_annotations'] = [x
    #                                   for x in variant['vep_annotations']
    #                                   if x['Gene'] == gene_id]
    #     add_consequence_to_variant(variant)
    #     variants.append(variant)
    # return variants
    pass


if __name__ == '__main__':
    db = scidbpy.connect()

    import pprint
    pp = pprint.PrettyPrinter(indent=2)
    pp.pprint(get_icd_significant(db, 'RH117'))
    pp.pprint(get_icd_info(db, 'RH117'))
    pp.pprint(get_variants_by_id(db, (1039381448,)))
    pp.pprint(get_variant(db, 1039381448))
    pp.pprint(get_variant_icd(db, 1039381448))
    pp.pprint(get_gene(db, 'ENSG00000101255'))
