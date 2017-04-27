GBE_FILES_DIRECTORY = '/home/scidb/GlobalBiobankEngine/gbe_data'
DBSNP_FILE = paste0(GBE_FILES_DIRECTORY, '/dbsnp142.txt.gz')
SITES_VCF_FILE = paste0(GBE_FILES_DIRECTORY, '/icd10ukbb.ukbiobank.merge.sort.vcf.gz')
ICD_STATS_FILES = system("find /home/scidb/GlobalBiobankEngine/gbe_data  -name 'c*.hybrid.gz' ", intern=T)
FILTER_FILE_1 = paste0(GBE_FILES_DIRECTORY, '/qc/UKBioBiLallfreqSNPexclude.dat')
FILTER_FILE_2 = paste0(GBE_FILES_DIRECTORY, '/qc/ukb_ukbl_low_concordance.dat')
ICD_INFO_FILE = paste0(GBE_FILES_DIRECTORY, '/icdstats/icdinfo.txt')
