"""
Utils for reading flat files that are loaded into database
"""
from __future__ import print_function
from __future__ import division
import re
import os
import traceback
from utils import *
import numpy

POPS = {
#    'AFR': 'African',
#    'AMR': 'Latino',
#    'EAS': 'East Asian',
    'FIN': 'European (Finnish)',
    'NFE': 'Non-Finnish European',
#    'SAS': 'South Asian',
    'OTH': 'European (Ashkenazi-Jewish)'
}


def get_base_coverage_from_file(base_coverage_file):
    """
    Read a base coverage file and return iter of dicts that look like:
    {
        'xpos': 1e9+1,
        'mean': 0.0,
        'median': 0.0,
        '1': 0.0,
        '5': 0.0,
        '10': 0.0,
        '15': 0.0,
        '20': 0.0,
        '25': 0.0,
        '30': 0.0,
        '50': 0.0,
        '100': 0.0,
    }
    """

    float_header_fields = ['odds_ratio', 'log_odds_ratio', 'pvalue', '-log10pvalue', 'flag', 'category']
    for line in base_coverage_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')
        d = {
            'xpos': get_xpos(fields[0], int(fields[1])),
            'pos': int(fields[1]),
        }
        for i, k in enumerate(float_header_fields):
            if i == 4:
                d[k] = fields[i+2]
            elif i == 5:
                if fields[i+2] == "NA":
                    continue
                if float_header_fields[i] == "category":
                    if csq_order_dict[fields[i+2]] <= csq_order_dict["frameshift_variant"]:
                        d[k] = 'lof_variant'
                    elif csq_order_dict[fields[i+2]] <= csq_order_dict["missense_variant"]:
                        d[k] = 'missense_variant'
                    elif csq_order_dict[fields[i+2]] <= csq_order_dict["synonymous_variant"]:
                        d[k] = 'synonymous_variant'
                    else:
                        d[k] = fields[i+2]
            else:
                d[k] = float(fields[i+2])
            # i+2
        yield d

def get_base_icdstats_from_file(base_icdstats_file):
    """
    Read a base coverage file and return iter of dicts that look like:
    {
        'xpos': 1e9+1,
        'pos' : value;
        'stats' : {'K50' : {'or' : value ; 'l95or' : value, 'u95or' : value, '-log10pvalue' : value, 'pvalue' : value}, 
        }
    }
    """

    float_header_fields = ['desc']
    for line in base_icdstats_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')
        d = {
              'xpos': get_xpos(fields[0], int(fields[1])),
            'pos': int(fields[1]),
             'stats' : []
        }
        # separate by ; 
        for icditem in fields[2].split(':'):
            d1 = {}
            for icdpair in icditem.split('/'):
                datkey = icdpair.split('=')[0]
                datvalue = icdpair.split('=')[1]
                d1[datkey] = datvalue
            d['stats'].append(d1)
            # i+2
        yield d

def get_icd_from_file(icdstats_file, filename):
    """
    Read an icd logistic file and return iter of dicts that look like:
    {
        'xpos': 1e9+1,
        'icd' : value,
        'pos' : value,
        'chrom' : value,
        'stats' : {'K50' : {'or' : value ; 'l95or' : value, 'u95or' : value, '-log10pvalue' : value, 'pvalue' : value}, 
        }
    }
    """
    fltdict = {}
    qcdat = open('../gbe_data/qc/ukb_ukbl_low_concordance.dat','r')
    for qcl in qcdat:
        qcl = qcl.rstrip()
        qcl = qcl.split()
        id = qcl[0]
        fltdict[id] = 'lowconcordance'
    qcdat.close()
    qcdat = open('../gbe_data/qc/UKBioBiLallfreqSNPexclude.dat','r')
    for qcl in qcdat:
        qcl = qcl.rstrip()
        qcl = qcl.split()
        id = qcl[0]
        fltdict[id] = 'exclude'
    qcdat.close()
    for line in icdstats_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split()
        if fields[6] != "ADD":
            continue
        if fields[2] in fltdict:
            continue
        if fields[11] == "NA":
            continue
        if float(fields[11]) == 0:
            continue
        if float(fields[9]) >= .5:
            continue
        fileidname = os.path.basename(filename[0]).split('.')
        if 'brainmri' in fileidname:
            prefix = 'BRMRI'
            intadd = 10
        elif 'additionalimaging' in fileidname:
            prefix = 'ADD'
            intadd = 20
        elif 'initialdata' in fileidname:
            prefix = 'INI'
            intadd = 30
        else:
            if len(fileidname[1].split('_FH2')) > 1:
                prefix = 'FH'
                intadd = 40
            elif len(fileidname[1].split('RH')) > 1:
                prefix = 'RH'
                intadd = 50
            elif len(fileidname[1].split('cancer')) > 1:
                prefix = 'cancer'
                intadd = 60
            elif len(fileidname[1].split('HC')) > 1:
                prefix = 'HC'
                intadd = 70
        d = {
              'icdind' : int(str(get_xpos(fields[0], int(fields[1]))) + str(intadd) + os.path.basename(filename[0]).split('.')[1].split('_')[0].strip()[1:].strip('RH').strip('FH').strip('cancer').strip('HC').split('_FH2')[0]),
              'icd' : prefix + os.path.basename(filename[0]).split('.')[1].strip('RH').strip('FH').strip('cancer').strip('HC').split('_FH2')[0],
                'xpos': get_xpos(fields[0], int(fields[1])),
              'affyid' : fields[2],
            'pos': int(fields[1]),
              'chrom' : fields[0],
             'stats' : []
        }
        d1 = {}
        d1['or'] = fields[8]
        d1['se'] = fields[9]
        d1['pvalue'] = float(fields[11])
        d1['lor'] = numpy.log(float(fields[8]))
        d1['log10pvalue'] = -float(numpy.log(float(fields[11]))/numpy.log(10))
        d1['l95or'] = numpy.exp(numpy.log(float(fields[8])) - 1.96*float(fields[9]))
        d1['u95or'] = numpy.exp(numpy.log(float(fields[8])) + 1.96*float(fields[9]))
        d['stats'].append(d1)
            # i+2
        yield d

def get_qt_from_file(qtstats_file, filename):
    """
    Read an icd logistic file and return iter of dicts that look like:
    {
        'xpos': 1e9+1,
        'icd' : value,
        'pos' : value,
        'chrom' : value,
        'stats' : {'K50' : {'or' : value ; 'l95or' : value, 'u95or' : value, '-log10pvalue' : value, 'pvalue' : value}, 
        }
    }
    """
    fltdict = {}
    qcdat = open('../gbe_data/qc/ukb_ukbl_low_concordance.dat','r')
    for qcl in qcdat:
        qcl = qcl.rstrip()
        qcl = qcl.split()
        id = qcl[0]
        fltdict[id] = 'lowconcordance'
    qcdat.close()
    qcdat = open('../gbe_data/qc/UKBioBiLallfreqSNPexclude.dat','r')
    for qcl in qcdat:
        qcl = qcl.rstrip()
        qcl = qcl.split()
        id = qcl[0]
        fltdict[id] = 'exclude'
    qcdat.close()
    for line in qtstats_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split()
        if fields[5] != "ADD":
            continue
        #CHROM POS ID REF ALT1 TEST OBS_CT BETA SE T_STAT P
        if fields[2] in fltdict:
            continue
        if fields[10] == "NA":
            continue
        se = float(fields[8])
        if se >= .5:
            continue
        icdn = os.path.basename(filename[0]).split('.')[1]
        icdcnt = {'1498' : 128300,
                  '1558' : 137558,
                  '20116' : 137558, 
                  '20159' :  32763,
                  '20240' : 30775,
                  '21001' :  137200}
        fileidname = os.path.basename(filename[0]).split('.')
        if 'brainmri' in fileidname:
            prefix = 'BRMRI'
            intadd = 10
        elif 'additionalimaging' in fileidname:
            prefix = 'ADD'
            intadd = 20
        elif 'initialdata' in fileidname:
            prefix = 'INI'
            intadd = 30
        else:
            if len(fileidname.split('_FH2')) > 1:
                prefix = 'FH'
                intadd = 40
            elif len(fileidname.split('RH')) > 1:
                prefix = 'RH'
                intadd = 50
            elif len(fileidname.split('cancer')) > 1:
                prefix = 'cancer'
                intadd = 60
            elif len(fileidname.split('HC')) > 1:
                prefix = 'HC'
                intadd = 70
        d = {
              'icdind' : int(str(get_xpos(fields[0], int(fields[1]))) + str(intadd) + os.path.basename(filename[0]).split('.')[1].split('_')[0].strip()[1:].strip('RH').strip('FH').strip('cancer').strip('HC').split('_FH2')[0]),
              'icd' : prefix + os.path.basename(filename[0]).split('.')[1].strip('RH').strip('FH').strip('cancer').strip('HC').split('_FH2')[0],
                'xpos': get_xpos(fields[0], int(fields[1])),
              'affyid' : fields[1],
            'pos': int(fields[1]),
              'chrom' : fields[0],
             'stats' : []
        }
        d1 = {}
        d1['or'] = float(fields[7])
        d1['se'] = float(fields[8])
        d1['pvalue'] = float(fields[10])
        d1['lor'] = float(fields[7])
        d1['log10pvalue'] = -float(numpy.log(d1['pvalue'])/numpy.log(10))
        d1['l95or'] = d1['lor'] - 1.96*d1['se']
        d1['u95or'] = d1['lor'] + 1.96*d1['se']
        d['stats'].append(d1)
            # i+2
        yield d

def get_icd_info_from_file(icdinfo_file):
    """
    Read an icd info file
        }
    }
    """
    icdf = open(icdinfo_file, 'r')
    for line in icdf:
        if line.startswith('#'):
            continue
        fields = line.rstrip().split('\t')
        name = fields[2]
        print(name)
        d = {
              'icd' : fields[0],
              'Name' : name.split('[')[0],
              'Case': fields[1]
        }
        yield d
        
def get_variants_from_sites_vcf(sites_vcf):
    """
    Parse sites VCF file and return iter of variant dicts
    sites_vcf is a file (gzipped), not file path
    """
    vep_field_names = None
    for line in sites_vcf:
        try:
            line = line.strip('\n')
#            if line.startswith('##INFO=<ID=CSQ'):
#                vep_field_names = line.split('Format: ')[-1].strip('">').split('|')
#                print(vep_field_names)
            if line.startswith('#'):
                continue
            # If we get here, it's a variant line
            if vep_field_names is None:
                vep_field_names = ['Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'Feature_type', 'Feature', 'BIOTYPE', 'EXON', 'INTRON', 'HGVSc', 'HGVSp', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation',  'DISTANCE', 'ALLELE_NUM','STRAND', 'VARIANT_CLASS', 'MINIMISED', 'SYMBOL_SOURCE', 'HGNC_ID', 'CANONICAL', 'TSL', 'APPRIS', 'CCDS', 'ENSP', 'SWISSPROT', 'TREMBL', 'UNIPARC',  'SIFT2','SIFT', 'PolyPhen', 'DOMAINS', 'HGVS_OFFSET', 'GMAF', 'AFR_MAF', 'AMR_MAF', 'EAS_MAF', 'EUR_MAF', 'SAS_MAF', 'AA_MAF', 'EA_MAF', 'ExAC_MAF', 'ExAC_Adj_MAF', 'ExAC_AFR_MAF', 'ExAC_AMR_MAF', 'ExAC_EAS_MAF', 'ExAC_FIN_MAF', 'ExAC_NFE_MAF', 'ExAC_OTH_MAF', 'ExAC_SAS_MAF', 'CLIN_SIG', 'SOMATIC', 'PHENO', 'PUBMED', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 'LoF', 'LoF_filter', 'LoF_flags', 'LoF_info']
            # This elegant parsing code below is copied from https://github.com/konradjk/loftee
            fields = line.split('\t')
            info_field = dict([(x.split('=', 1)) if '=' in x else (x, x) for x in re.split(';(?=\w)', fields[7])])
            consequence_array = info_field['CSQ'].split(',') if 'CSQ' in info_field else []
            annotations = [dict(zip(vep_field_names, x.split('|'))) for x in consequence_array]
            coding_annotations = [ann for ann in annotations if 'Feature' in ann and ann['Feature'].startswith('ENST')]
            alt_alleles = fields[4].split(',')
            # different variant for each alt allele
            for i, alt_allele in enumerate(alt_alleles):
                vep_annotations = [ann for ann in coding_annotations]
#                print([ann['DISTANCE'] for ann in coding_annotations])
 #               print([ann['CANONICAL'] for ann in coding_annotations])
                #exit()
                # Variant is just a dict
                # Make a copy of the info_field dict - so all the original data remains
                # Add some new keys that are allele-specific
                pos, ref, alt = get_minimal_representation(fields[1], fields[3], alt_allele)
                variant = {}
                variant['chrom'] = fields[0]
                variant['pos'] = pos
                variant['rsid'] = fields[2]
                variant['xpos'] = get_xpos(variant['chrom'], variant['pos'])
                variant['ref'] = ref
                variant['alt'] = alt
                variant['xstart'] = variant['xpos']
                variant['xstop'] = variant['xpos'] + len(variant['alt']) - len(variant['ref'])
                variant['variant_id'] = '{}-{}-{}-{}'.format(variant['chrom'], variant['pos'], variant['ref'], variant['alt'])
                variant['orig_alt_alleles'] = [
                    '{}-{}-{}-{}'.format(variant['chrom'], *get_minimal_representation(fields[1], fields[3], x))
                    for x in alt_alleles
                ]
                variant['site_quality'] = fields[5]
                variant['filter'] = fields[6]
                variant['vep_annotations'] = vep_annotations

                variant['minicd'] = info_field['minicd']
                variant['minpval'] = float(info_field['minpval'])
                variant['minor'] = float(info_field['minor'])
                variant['minl10pval'] = float(info_field['minl10pval'])
#                if not variant['allele_count'] and variant['filter'] == 'PASS': variant['filter'] = 'AC_Adj0' # Temporary filter
                if 'EXAC_NFE' not in info_field:
                    variant['exac_nfe'] = str('-9')
                    variant['allele_freq'] = str('-9')
                    variant['allele_num'] = str('-9')
                else:
                    variant['exac_nfe'] = info_field['EXAC_NFE'].split(',')[0]
                    variant['allele_freq'] = info_field['EXAC_NFE'].split(',')[0]
                    variant['allele_num'] = str('-9')
                variant['genes'] = list({annotation['Gene'] for annotation in vep_annotations})
                variant['transcripts'] = list({annotation['Feature'] for annotation in vep_annotations})                   
                yield variant
        except Exception:
            print("Error parsing vcf line: " + line)
            traceback.print_exc()
            break


def get_canonical_transcripts(canonical_transcript_file):
    for line in canonical_transcript_file:
        gene, transcript = line.strip().split()
        yield gene, transcript


def get_omim_associations(omim_file):
    for line in omim_file:
        fields = line.strip().split('\t')
        if len(fields) == 4:
            yield fields
        else:
            yield None


def get_genes_from_gencode_gtf(gtf_file):
    """
    Parse gencode GTF file;
    Returns iter of gene dicts
    """
    for line in gtf_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')

        if fields[2] != 'gene':
            continue

        chrom = fields[0][3:]
        start = int(fields[3]) + 1  # bed files are 0-indexed
        stop = int(fields[4]) + 1
        info = dict(x.strip().split() for x in fields[8].split(';') if x != '')
        info = {k: v.strip('"') for k, v in info.items()}
        gene_id = info['gene_id'].split('.')[0]

        gene = {
            'gene_id': gene_id,
            'gene_name': info['gene_name'],
            'gene_name_upper': info['gene_name'].upper(),
            'chrom': chrom,
            'start': start,
            'stop': stop,
            'strand': fields[6],
            'xstart': get_xpos(chrom, start),
            'xstop': get_xpos(chrom, stop),
        }
        yield gene


def get_transcripts_from_gencode_gtf(gtf_file):
    """
    Parse gencode GTF file;
    Returns iter of transcript dicts
    """
    for line in gtf_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')

        if fields[2] != 'transcript':
            continue

        chrom = fields[0][3:]
        start = int(fields[3]) + 1  # bed files are 0-indexed
        stop = int(fields[4]) + 1
        info = dict(x.strip().split() for x in fields[8].split(';') if x != '')
        info = {k: v.strip('"') for k, v in info.items()}
        transcript_id = info['transcript_id'].split('.')[0]
        gene_id = info['gene_id'].split('.')[0]

        gene = {
            'transcript_id': transcript_id,
            'gene_id': gene_id,
            'chrom': chrom,
            'start': start,
            'stop': stop,
            'strand': fields[6],
            'xstart': get_xpos(chrom, start),
            'xstop': get_xpos(chrom, stop),
        }
        yield gene


def get_exons_from_gencode_gtf(gtf_file):
    """
    Parse gencode GTF file;
    Returns iter of transcript dicts
    """
    for line in gtf_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')

        if fields[2] not in ['exon', 'CDS', 'UTR']:
            continue

        chrom = fields[0][3:]
        feature_type = fields[2]
        start = int(fields[3]) + 1  # bed files are 0-indexed
        stop = int(fields[4]) + 1
        info = dict(x.strip().split() for x in fields[8].split(';') if x != '')
        info = {k: v.strip('"') for k, v in info.items()}
        transcript_id = info['transcript_id'].split('.')[0]
        gene_id = info['gene_id'].split('.')[0]

        exon = {
            'feature_type': feature_type,
            'transcript_id': transcript_id,
            'gene_id': gene_id,
            'chrom': chrom,
            'start': start,
            'stop': stop,
            'strand': fields[6],
            'xstart': get_xpos(chrom, start),
            'xstop': get_xpos(chrom, stop),
        }
        yield exon


def get_dbnsfp_info(dbnsfp_file):
    """
    Parse dbNSFP_gene file;
    Returns iter of transcript dicts
    """
    header = dbnsfp_file.next().split('\t')
    fields = dict(zip(header, range(len(header))))
    for line in dbnsfp_file:
        line = line.split('\t')
        other_names = line[fields["Gene_old_names"]].split(';') if line[fields["Gene_old_names"]] != '.' else []
        if line[fields["Gene_other_names"]] != '.':
            other_names.extend(line[fields["Gene_other_names"]].split(';'))
        gene_info = {
            'gene_name': line[fields["Gene_name"]],
            'ensembl_gene': line[fields["Ensembl_gene"]],
            'gene_full_name': line[fields["Gene_full_name"]],
            'gene_other_names': other_names
        } 
        yield gene_info


def get_snp_from_dbsnp_file(dbsnp_file):
    for line in dbsnp_file:
        fields = line.split('\t')
        rsid = int(fields[0])
        chrom = fields[1].rstrip('T')
        if chrom == 'PAR': continue
        start = int(fields[2]) + 1
        snp = {
            'xpos': get_xpos(chrom, start),
            'rsid': rsid
        }
        yield snp
