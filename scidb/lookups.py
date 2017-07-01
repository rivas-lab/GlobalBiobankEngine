import exceptions
import scidbpy


import config


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
                   icd.icd_id,
                   icd_index.icd_id),
        pvalue < 0.01);
    """
    return numpy2dict(
        db.iquery(
            config.ICD_LOOKUP_QUERY.format(icd_id=icd_id, cutoff=cutoff),
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
                 icd_info.icd_id,
                 icd_index.icd_id);
    """
    return numpy2dict(
        db.iquery(
            config.ICD_INFO_LOOKUP_QUERY.format(icd_id=icd_id),
            schema=config.ICD_INFO_LOOKUP_SCHEMA_INST,
            fetch=True))


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
    print(variant_ids)
    if len(variant_ids) > 1:
        raise exceptions.NotImplementedError()
    return numpy2dict(
        db.iquery(
            config.VARIANT_LOOKUP_QUERY.format(xpos=variant_ids[0]),
            schema=config.VARIANT_LOOKUP_SCHEMA_INST,
            fetch=True))


if __name__ == '__main__':
    db = scidbpy.connect()

    import pprint
    pp = pprint.PrettyPrinter(indent=2)
    pp.pprint(get_icd_significant(db, 'RH117')[:1])
    pp.pprint(get_icd_info(db, 'RH117'))
    pp.pprint(get_variants_by_id(db, (1039381448,))[:1])
