import os.path

GBE_DATA_PATH = '/home/scidb/GlobalBiobankEngine/gbe_data'

QC_FILTER_PATH = os.path.join(GBE_DATA_PATH, 'qc')
QC_FILTER_FILES = (
    {'filename': os.path.join(QC_FILTER_PATH, 'UKBioBiLallfreqSNPexclude.dat'),
     'header': 1},
    {'filename': os.path.join(QC_FILTER_PATH, 'ukb_ukbl_low_concordance.dat')}
)
