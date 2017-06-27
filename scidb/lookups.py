import scidbpy


import config


def get_icd_significant(db, icd_id, cutoff=0.01):
    """
    e.g.,
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


if __name__ == '__main__':
    db = scidbpy.connect()
    print(get_icd_significant(db, 'RH117'))
