# SciDB BioBank API Quick Start

Activate Python virtual environment and start IPython:

    ~> source /opt/biobank-api/bin/activate
    (biobank-api)~> ipython

Disable warnings (optional):

    import requests
    requests.packages.urllib3.disable_warnings(
        requests.packages.urllib3.exceptions.SNIMissingWarning)
    requests.packages.urllib3.disable_warnings(
        requests.packages.urllib3.exceptions.InsecurePlatformWarning)

    import warnings
    warnings.filterwarnings('ignore', category = UserWarning)
    warnings.filterwarnings('ignore', category = FutureWarning)

Import BioBank API and connect to SciDB:

    from scidbbiobank import connect
    bb = connect(scidb_auth=('scidbadmin', 'Paradigm4'), namespace='RIVAS_HG19')

Get genes matching `*RP11-150*`:

    bb.get_genes(gene_name='RP11-150', exact_match=False)

Get genes overlapping with specified region:

    bb.get_genes_in_region(chromosome=2, start=202201300, end=202201800)

Get transcripts and exons using gene or transcripts ensembl IDs:

    bb.get_transcripts(namespace='RIVAS_HG38', gene_eid='ENSG00000071626')

    bb.get_transcripts(transcript_eid=['ENST00000233078', 'ENST00000250863'])

    bb.get_exons(gene_eid='ENSG00000071626', transcript_eid='ENST00000233078')

Get variant fields:

    bb.get_variant_fields()

Get variants using various criteria:

    bb.get_variants(chromosome=22, position=(32334104, 38877475))

    bb.get_variants(chromosome=22, start=32300000, end=32400000)

    bb.get_variants(chromosome=(21, 22), start=32300000, end=32400000)

    bb.get_variants(rsid=['rs34221567', 'rs148534902', 'rs200712517'])

    bb.get_variants(chromosome=22,
        consequence_contains='frameshift', impact='SINGLE_EXON')

    bb.get_variants(chromosome=22,
        consequence_contains='frameshift', miss_range=(0.8, None))

    bb.get_variants(gene_name=['RP11-150L8.3', 'RP11-150D20.5'], pad=10000)

List association sets:

    bb.list_association_sets()
    assocset = str(bb.list_association_sets()['name'][0])

Get association set data:

    bb.get_association_data(association_set=assocset,
        chromosome=22, start=32300000, end=32400000, pvalue_max=.0000001)

    bb.get_association_data(association_set=assocset,
        gene_name='RP11-150D20.5', pad=10000, pvalue_max=.000001)

    bb.get_association_data(association_set=assocset,
        gene_name='RP11-150D20.5', pad=10000, pvalue_max=.000001, field_id = (0))

Any criteria from `get_variants` can be used in `get_association_data`
as well.

Get phenotype fields:

    df = bb.get_phenotype_fields(association_set=assocset)
    df[df['description'] == 'asthma']
        instance_id  value_no  field_id  title value_type description                          notes
        1719           21        12       499  HC382       None      asthma  short_name:Asthma;category:HC
    str(list(df[df['description'] == 'asthma']['title'])[0])
        'HC382'
    int(df[df['description'] == 'asthma']['field_id'])
        499
    bb.get_phenotype_fields(association_set='ARRAY', name_kw='BIN1')
    bb.get_phenotype_fields(association_set='ARRAY', name_kw=('BIN12', 'BIN21'))

Delete association fields:

    bb.detele_association_data(association_set='FOO',
                               field_id=100)

Delete entire association set *use with caution*:

    bb.delete_association_set(association_set='FOO')

Leave virtual environment after exiting IPython:

    (biobank-api)~> deactivate
