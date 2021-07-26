# BioBank Streaming Example

    variant_ids = {
        0 : bb.get_variants(gene_name='SYNC')['variant_id'].values,
        1 : bb.get_variants(gene_name='RP11-114B7.6')['variant_id'].values}

    field_id = bb.get_phenotype_fields(
        association_set='UK_Biobank_Array_European',
        name_kw=('RH24', 'RH25'))['field_id'].values

    bb.stream(
        'RIVAS_HG19', 'UK_Biobank_Array_European', variant_ids, field_id)
       field_id  list_id   key value title
    0         0        1  noop  None  RH24
    1         0        0  noop  None  RH24
    2         1        1  noop  None  RH25
    3         1        0  noop  None  RH25


    def simple(df):
        return {
            'min_odds_ratio': min(df['odds_ratio']),
            'avg_se': sum(df['se']) / len(df)}

    bb.stream(
        'RIVAS_HG19', 'UK_Biobank_Array_European', variant_ids, field_id, simple)
       field_id  list_id             key           value title
    0         0        1  min_odds_ratio        0.909121  RH24
    1         0        1          avg_se       0.0578869  RH24
    2         0        0  min_odds_ratio        0.699871  RH24
    3         0        0          avg_se    0.3817177125  RH24
    4         1        1  min_odds_ratio         0.69296  RH25
    5         1        1          avg_se  0.634265333333  RH25
    6         1        0  min_odds_ratio       0.0236586  RH25
    7         1        0          avg_se    0.9201491875  RH25
