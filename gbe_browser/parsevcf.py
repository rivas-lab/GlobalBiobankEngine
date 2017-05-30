from __future__ import print_function
from __future__ import division

import os
import sys
import numpy as np
import gzip 
import glob
from parsing import *
from utils import *



chroffdict = {}

chrmap = {}
chrst = {}
resd = {}

affydict = {}

for file in glob.glob('../gbe_data/icdassoc/hybrid/*hybrid.rewrite.gz'):
    fr = gzip.open(file,'rb')
    print(file)
    fileidname = os.path.basename(file).split('.')
    if 'brainmri' in fileidname:
        prefix = 'BRMRI'
    elif 'additionalimaging' in fileidname:
        prefix = 'ADD'
    elif 'initialdata' in fileidname:
        prefix = 'INI'
    else:
        if len(fileidname[1].split('_FH2')) > 1:
            prefix = 'FH'
        elif len(fileidname[1].split('RH')) > 1:
            prefix = 'RH'
        elif len(fileidname[1].split('cancer')) > 1:
            prefix = 'cancer'
        elif len(fileidname[1].split('HC')) > 1:
            prefix = 'HC'
    icd = prefix + os.path.basename(file).split('.')[1].strip('RH').strip('FH').strip('cancer').strip('HC').split('_FH2')[0]
    for line in fr:
        if line[0] == "#":
            continue
        line = line.rstrip()
        line = line.split()
        chrom = line[0]
        pos = line[1] 
        affyid = line[2]
        if line[6] != "ADD":
            continue
        if line[11] == "NA":
            continue
        fields = line
        if fields[11] == "NA":
            continue
        if float(fields[9]) > 1:
            continue
        if float(fields[11]) == 0:
            continue
        pval = float(line[11])
        l10pval = -np.log(pval)/np.log(10)
        oddsratio = float(line[8])
        ser = float(line[9])
        if ser >= .5:
            continue
        l95or = np.exp(np.log(oddsratio) - 1.96*ser) 
        l95or  = '{0:.2f}'.format(float(l95or))
        u95or = np.exp(np.log(oddsratio) + 1.96*ser)
        u95or  = '{0:.2f}'.format(float(u95or))
        oddsratio  = '{0:.2f}'.format(float(oddsratio))
        l10pval  = '{0:.2f}'.format(float(l10pval))
        if (affyid,'pval') in chrst:
            if pval <= .01:
                chrst[affyid,'pval'].append(pval) 
                chrst[affyid,'oddsr'].append(oddsratio) 
                chrst[affyid,'l10pval'].append(l10pval) 
                chrst[affyid,'icd'].append(icd) 
        else:
            chrst[affyid,'pval'] = [pval] 
            chrst[affyid,'oddsr'] = [oddsratio] 
            chrst[affyid,'l10pval'] = [l10pval] 
            chrst[affyid,'icd'] = [icd]  
    fr.close()
        

for file in glob.glob('../gbe_data/icdassoc/hybrid/*linear.rewrite.gz'):
    fr = gzip.open(file,'rb')
    fileidname = os.path.basename(file).split('.')
    if 'brainmri' in fileidname:
        prefix = 'BRMRI'
    elif 'additionalimaging' in fileidname:
        prefix = 'ADD'
    elif 'initialdata' in fileidname:
        prefix = 'INI'
    else:
        if len(fileidname[1].split('_FH2')) > 1:
            prefix = 'FH'
        elif len(fileidname[1].split('RH')) > 1:
            prefix = 'RH'
        elif len(fileidname[1].split('cancer')) > 1:
            prefix = 'cancer'
        elif len(fileidname[1].split('HC')) > 1:
            prefix = 'HC'
    icd = prefix + os.path.basename(file).split('.')[1].strip('RH').strip('FH').strip('cancer').strip('HC').split('_FH2')[0]
    for line in fr:
        if line[0] == "#":
            continue
        line = line.rstrip()
        line = line.split()
        chrom = line[0]
        pos = line[1] 
        affyid = line[2]
        if line[5] != "ADD":
            continue
        if line[7] == "NA":
            continue
        beta = float(line[7])
        ser = float(line[8])
        if ser >= .5:
            continue
        l95beta = beta - 1.96*ser
        u95beta = beta + 1.96*ser 
        u95beta = '{0:.2f}'.format(float(u95beta))
        l95beta = '{0:.2f}'.format(float(l95beta))
        beta  = '{0:.2f}'.format(float(beta))
        pval = float(line[10])
        l10pval = -np.log(pval)/np.log(10)
        oddsratio = beta
        l95or = l95beta
        u95or = u95beta
        if (affyid,'pval') in chrst:
            if pval <= .01:
                chrst[affyid,'pval'].append(pval) 
                chrst[affyid,'oddsr'].append(oddsratio) 
                chrst[affyid,'l10pval'].append(l10pval) 
                chrst[affyid,'icd'].append(icd) 
        else:
            chrst[affyid,'pval'] = [pval] 
            chrst[affyid,'oddsr'] = [oddsratio] 
            chrst[affyid,'l10pval'] = [l10pval] 
            chrst[affyid,'icd'] = [icd]  
    fr.close()
        

exacnfe = {}
seendict = {}

vcfr = gzip.open(sys.argv[2],'rb')

vcfout = gzip.open('icd10ukbb.' + sys.argv[2],'wb')
covout = gzip.open('../gbe_data/icdstats/Panel2016.all.coverage.txt.gz','wb')

for line in vcfr:
    line = line.rstrip('\n')
    if line.startswith('##INFO=<ID=CSQ'):
        vep_field_names = line.split('Format: ')[-1].strip('">').split('|')
    if line[0] == "#":
        pass
    else:
        line = line.rstrip()
        line = line.split('\t')
        affyid = line[6]        
        chroffset = str(line[0]) + ":" + str(line[1])
        seendict[chroffset] = chroffset 

vcfr.close()

vcff = gzip.open(sys.argv[1],'rb')
for line in vcff:
    line = line.strip('\n')
    if line.startswith('##INFO=<ID=CSQ'):
        vep_field_names = line.split('Format: ')[-1].strip('">').split('|')
    if line[0] == "#":
        pass
    else:
        line = line.rstrip()
        line = line.split('\t')
        chroffset = line[0] + ":" + line[1] 
        if chroffset not in seendict:
            continue
        info_field = dict([(x.split('=', 1)) if '=' in x else (x, x) for x in re.split(';(?=\w)', line[7])])
        if 'AF_NFE' not in info_field:
            exacnfe[chroffset] = str('-9')
        else:
            exacnfe[chroffset] = info_field['AF_NFE'].split(',')[0]
vcff.close()

vcfr = gzip.open(sys.argv[2],'rb')

vcfout = gzip.open('icd10ukbb.' + sys.argv[2],'wb')
covout = gzip.open('../gbe_data/icdstats/Panel2016.all.coverage.txt.gz','wb')
seen = {}
for line in vcfr:
    line = line.rstrip('\n')
    if line.startswith('##INFO=<ID=CSQ'):
        vep_field_names = line.split('Format: ')[-1].strip('">').split('|')
        print(line,file = vcfout)
    if line[0] == "#":
        pass
    else:
        line = line.rstrip()
        line = line.split('\t')
        affyid = line[6]        
        chroffset = str(line[0]) + ":" + str(line[1])
        if chroffset in seen:
            continue
        else:
            seen[chroffset] = chroffset
        if chroffset in exacnfe:
            exacfreq = exacnfe[chroffset]
        else:
            exacfreq = '-9'
        if (affyid,'pval') in chrst:
            pvalarr = [float(chrst[affyid,'pval'][i]) for i in range(0,len(chrst[affyid,'pval']))]
            ind = pvalarr.index(min(pvalarr))
            minpval = pvalarr[ind]
            minicd = chrst[affyid,'icd'][ind]
            minor = chrst[affyid,'oddsr'][ind]
            minl10pval = chrst[affyid,'l10pval'][ind]
            info_field = dict([(x.split('=', 1)) if '=' in x else (x, x) for x in re.split(';(?=\w)', line[7])])
            info = 'EXAC_NFE=' + str(exacfreq) + ';' + 'minicd=' + str(minicd) + ';' + 'minpval=' + str(minpval) + ';' + 'minor=' + str(minor) + ';' + 'minl10pval=' + str(minl10pval) + ';' + line[7]
            consequence_array = info_field['CSQ'].split(',') if 'CSQ' in info_field else []
            annotations = [dict(zip(vep_field_names, x.split('|'))) for x in consequence_array if len(vep_field_names) == len(x.split('|'))]
            coding_annotations = [ann for ann in annotations if ann['Feature'].startswith('ENST')]
            csq = order_vep_by_csq(coding_annotations)
            worst_csq = worst_csq_with_vep(csq)
            if worst_csq is None:
                worst = 'intergenic'
            else:
                worst = worst_csq['major_consequence']
            print(line[0],line[1],line[2],line[3],line[4],'.','PASS',info,sep = '\t', file = vcfout)
            print(line[0],line[1],minor,str(np.log(float(minor))),minpval,minl10pval,'PASS',worst,file = covout,sep = '\t')



vcfr.close()


covout.close()
vcfout.close()

