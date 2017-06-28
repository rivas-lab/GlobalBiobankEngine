import scidbpy


import config


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
    return db.iquery(
        config.ICD_LOOKUP_QUERY.format(icd_id=icd_id, cutoff=cutoff),
        schema=config.ICD_LOOKUP_SCHEMA,
        fetch=True,
        as_dataframe=True)


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
    return db.iquery(
        config.ICD_INFO_LOOKUP_QUERY.format(icd_id=icd_id),
        schema=config.ICD_INFO_LOOKUP_SCHEMA,
        fetch=True,
        as_dataframe=True)


if __name__ == '__main__':
    db = scidbpy.connect()
    print(get_icd_significant(db, 'RH117'))
    print(get_icd_info(db, 'RH117'))
