from __future__ import division
from __future__ import print_function
import dash
import dash_core_components as dcc
import dash_html_components as html
import sys
import itertools
import collections
import io
import json
import os
import re
import gzip
from parsing import *
import lookups
import random
import numpy
import numpy as np
import numpy.matlib
import pandas
from utils import *
import utils
from flask import Flask, request, session, g, redirect, url_for, abort, render_template, flash, jsonify
from flask_compress import Compress
from flask_oauth import OAuth
from flask_errormail import mail_on_500
from flask import Response
from collections import defaultdict
from werkzeug.contrib.cache import SimpleCache
from multiprocessing import Process
import glob
import traceback
import time
import csv
import logging
from scipy import stats
from flask_wtf import Form, RecaptchaField
from wtforms import StringField, BooleanField
from wtforms.validators import DataRequired
### New SciDB requirements
import requests
requests.packages.urllib3.disable_warnings(
    requests.packages.urllib3.exceptions.SNIMissingWarning)
requests.packages.urllib3.disable_warnings(
    requests.packages.urllib3.exceptions.InsecurePlatformWarning)

import warnings

from scidbbiobank import connect
###### End SciDB requirements


from dash_apps import *
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go


ADMINISTRATORS = (
    'gbe.help@gmail.com',
)

app = Flask(__name__)
app.debug=True

# Set random seed

mail_on_500(app, ADMINISTRATORS)
COMPRESS_MIMETYPES = ['text/html', 'text/css', 'application/json']
COMPRESS_LEVEL = 6
COMPRESS_MIN_SIZE = 500
Compress(app)
app.config['COMPRESS_DEBUG'] = True

cache = SimpleCache()


GBE_FILES_DIRECTORY = '../gbe_data/'
REGION_LIMIT = 1E5
EXON_PADDING = 50
# Load default config and override config from an environment variable
#test
app.config.update(dict(
    ## Set SCIDB_URL='http://localhost:8080' environment variable
    # SCIDB_URL='http://localhost:8080',
    DEBUG=False,
    SECRET_KEY='development key',
    LOAD_DB_PARALLEL_PROCESSES = 8,  # contigs assigned to threads, so good to make this a factor of 24 (eg. 2,3,4,6,8)
    SITES_VCFS=glob.glob(os.path.join(os.path.dirname(__file__), GBE_FILES_DIRECTORY, 'icd10ukbb.ukbiobank.combined.vcf.gz')),
    GENCODE_GTF=os.path.join(os.path.dirname(__file__), GBE_FILES_DIRECTORY, 'gencode.gtf.gz'),
    CANONICAL_TRANSCRIPT_FILE=os.path.join(os.path.dirname(__file__), GBE_FILES_DIRECTORY, 'canonical_transcripts.txt.gz'),
    OMIM_FILE=os.path.join(os.path.dirname(__file__), GBE_FILES_DIRECTORY, 'omim_info.txt.gz'),
    BASE_COVERAGE_FILES=glob.glob(os.path.join(os.path.dirname(__file__), GBE_FILES_DIRECTORY, 'coverage', 'Panel2016.*.coverage.txt.gz')),
    BASE_ICDSTATS_FILES=glob.glob(os.path.join(os.path.dirname(__file__), GBE_FILES_DIRECTORY, 'icdstats', 'Panel*.icdstats.txt.gz')),
    ICD_INFO_FILE=os.path.join(os.path.dirname(__file__), GBE_FILES_DIRECTORY, 'icdstats', 'icdinfo.txt'),
    ICD_STATS_FILES=glob.glob(os.path.join(os.path.dirname(__file__), GBE_FILES_DIRECTORY, 'icdassoc','hybrid','*c*.hybrid.rewritewna.gz')),
    QT_STATS_FILES=glob.glob(os.path.join(os.path.dirname(__file__), GBE_FILES_DIRECTORY, 'icdassoc','hybrid','*c*.linear.rewritewna.gz')),
    DBNSFP_FILE=os.path.join(os.path.dirname(__file__), GBE_FILES_DIRECTORY, 'dbNSFP2.6_gene.gz'),
    # This is not supported in GBE, We have dbsnp150
    # How to get a snp141.txt.bgz file:
    #   wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp141.txt.gz
    #   zcat snp141.txt.gz | cut -f 1-5 | bgzip -c > snp141.txt.bgz
    #   tabix -0 -s 2 -b 3 -e 4 snp141.txt.bgz

    #   wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/database/organism_data/b142_SNPChrPosOnRef_105.bcp.gz
    #   zcat b142_SNPChrPosOnRef_105.bcp.gz | awk '$3 != ""' | perl -pi -e 's/ +/\t/g' | sort -k2,2 -k3,3n | bgzip -c > dbsnp142.txt.bgz
    #   tabix -s 2 -b 3 -e 3 dbsnp142.txt.bgz
    DBSNP_FILE=os.path.join(os.path.dirname(__file__), GBE_FILES_DIRECTORY, 'dbsnp150.txt.gz')
))
GENE_CACHE_DIR = os.path.join(os.path.dirname(__file__), 'gene_cache')
GENES_TO_CACHE = {l.strip('\n') for l in open(os.path.join(os.path.dirname(__file__), 'genes_to_cache.txt'))}
### New SciDB connect
f = open('/biobankengine/app/static/.credentials','r').readlines()
uname = f[0].rstrip()
pwd = f[1].rstrip()
DB = connect(scidb_url=os.getenv('SCIDB_URL',None), scidb_auth=(uname, pwd), namespace='RIVAS_HG19')
### End SciDB connect



def create_cache():
    """
    This is essentially a compile step that generates all cached resources.
    Creates files like autocomplete_entries.txt
    Should be run on every redeploy.
    """
    # create autocomplete_entries.txt
    print >> sys.stderr, "Done!"



def get_db(name_space):
    """Opens a new database connection if there is none yet for the
    current application context.

    """
    # if not hasattr(g, 'db_conn'):
    #     g.db_conn = connect_db()
    # return g.db_conn
    DB = connect(scidb_url=os.getenv('SCIDB_URL',None), scidb_auth=(uname, pwd), namespace=name_space)
    DB.set_limit(15000)
    return DB

@app.route('/', methods = ['GET', 'POST'])
def homepage():
    if request.method == 'POST':
        namespace = request.form['functionassocset']
        db = get_db(namespace)
        nsdesc = db.list_association_sets()
        print(db.list_association_sets())
        ns = [{'id': namespace, 'name' : ' '.join(nsdesc['name'][0].split('_'))}]
        return render_template(
            'search_page_ns.html',
            ns=ns,
            namespace=namespace,
            formset=namespace
            )
    else:
        namespace = 'RIVAS_HG19'
        db = get_db(namespace)
        nsdesc = db.list_association_sets()
        print(db.list_association_sets())
        ns = [{'id': namespace, 'name' : ' '.join(nsdesc['name'][0].split('_'))}]
        return render_template(
            'search_page_ns.html',
            ns=ns,
            namespace=namespace,
            formset=namespace
            )
    return render_template('search_page.html')

@app.route('/search', methods = ['GET', 'POST'])
def search_page():
    if request.method == 'POST':
        namespace = request.form['functionassocset']
        db = get_db(namespace)
        nsdesc = db.list_association_sets()
        print(db.list_association_sets())
        ns = [{'id': namespace, 'name' : ' '.join(nsdesc['name'][0].split('_'))}]
        return render_template(
            'search_page_ns.html',
            ns=ns,
            namespace=namespace,
            formset=namespace
            )
    return render_template('search_page.html')

@app.route('/search/<namespace>', methods = ['GET', 'POST'])
def searchnamespace_page(namespace):
#    if not check_credentials():
#        return redirect(url_for('login'))
    if request.method == 'POST':
        #  print(dir(request))
        #  print(request.form.keys())
        namespace = request.form['functionassocset']
        db = get_db(namespace)
        nsdesc = db.list_association_sets()
        print(db.list_association_sets())
        ns = [{'id': namespace, 'name' : ' '.join(nsdesc['name'][0].split('_'))}]
        return render_template(
            'search_page_ns.html',
            ns=ns,
            namespace=namespace,
            formset=namespace
            )
    else:
        db = get_db(namespace)
        nsdesc = db.list_association_sets()
        print(db.list_association_sets())
        ns = [{'id': namespace, 'name' : ' '.join(nsdesc['name'][0].split('_'))}]
        return render_template(
            'search_page_ns.html',
            ns=ns,
            namespace=namespace,
            formset=namespace
            )
    return render_template('search_page.html')


@app.route('/<namespace>/autocomplete/<query>')
def awesome_autocomplete(namespace, query):
    if not hasattr(g, 'autocomplete_strings'):
        g.autocomplete_strings = [s.strip() for s in open(os.path.join(os.path.dirname(__file__), namespace + '_autocomplete_strings.txt'))]
    suggestions = lookups.get_awesomebar_suggestions(g, query)
    return Response(json.dumps([{'value': s} for s in suggestions]),  mimetype='application/json')


@app.route('/<namespace>/awesome')
def awesome(namespace):
    db = get_db(namespace)
    query = request.args.get('query')
    datatype, identifier = lookups.get_awesomebar_result(db, query)
    print("Searched for %s: %s" % (datatype, identifier))
    if datatype == 'gene':
        return redirect('/{0}/gene/{1}'.format(namespace,identifier))
    elif datatype == 'transcript':
        return redirect('/{0}/transcript/{1}'.format(namespace,identifier))
    elif datatype == 'dbsnp_variant_set':
        return redirect('/{0}/dbsnp/{1}'.format(namespace,identifier))
    elif datatype == 'variant':
        return redirect('/{0}/variant/{1}'.format(namespace,identifier))
    elif datatype == 'gbe':
        return redirect('/{0}/coding/{1}'.format(namespace,identifier))
    elif datatype == 'region':
        return redirect('/{0}/region/{1}'.format(namespace,identifier))
    elif datatype == 'error':
        return redirect('/{0}/error/{1}'.format(namespace,identifier))
    elif datatype == 'not_found':
        return redirect('/{0}/not_found/{1}'.format(namespace,identifier))
    else:
        raise Exception


@app.route('/publications')
def publications_page():
    namespace = 'RIVAS_HG19'
    try:
        return render_template('publications.html',
                               namespace = namespace)
    except Exception as e:
        print('Unknown Error=', traceback.format_exc())
        abort(404)


@app.route('/<namespace>/variant/<variant_str>')
def variant_icd_page(namespace, variant_str):
    db = get_db(namespace)
    try:
        chrom, pos, ref, alt = variant_str.split('-')
        chrom = int(chrom)
        pos = int(pos)
        refl = [ref]
        altl = [alt]
        # pos, ref, alt = get_minimal_representation(pos, ref, alt)
        variant = lookups.get_variant_ann_by_chrom_pos(db, chrom, pos, ref, alt)
        if variant is None:
            variant = {
                'chrom': chrom,
                'pos': pos,
                'ref' : ref,
                'alt' : alt,
                'xpos': get_xpos(chrom, int(pos))
            }
        consequences = None
        ordered_csqs = None
        #chrom       pos ref alt variant_identity major_consequence          category        gene_name gene_symbol                         HGVSp                       HGVSc       consequence  all_filters                                        annotations       maf  ukbb_freq    ld        rsid
        variant = variant.to_dict('records')[0]
        if 'annotations' in variant:
            vep_annotation = str(variant['annotations'])
            info_field = dict([(x.split('=', 1)) if '=' in x else (x, x) for x in re.split(';(?=\w)', vep_annotation)])
            consequence_array = info_field['CSQ'].split(',') if 'CSQ' in info_field else []
            annotations = [dict(zip(utils.vep_field_names, x.split('|'))) for x in consequence_array]
            coding_annotations = [ann for ann in annotations if 'Feature' in ann and ann['Feature'].startswith('ENST')]
            vep_annotations = [ann for ann in coding_annotations]
            variant['vep_annotations'] = vep_annotations
            variant['filter'] =  'PASS'
            #if variant['all_filters'] != 1 and variant['all_filters'] != 2 and variant['all_filters'] != 3 else 'FAIL'
            # variant['vep_annotations'] = order_vep_by_csq(variant['vep_annotations'])  # Adds major_consequence
           # ordered_csqs = [x['major_consequence'] for x in variant['vep_annotations']]
           # ordered_csqs = reduce(lambda x, y: ','.join([x, y]) if y not in x else x, ordered_csqs, '').split(',') # Close but not quite there
           # consequences = defaultdict(lambda: defaultdict(list))
           # for annotation in variant['vep_annotations']:
           #     annotation['HGVS'] = get_proper_hgvs(annotation)
           #     consequences[annotation['major_consequence']][annotation['Gene']].append(annotation)
            variant['vep_annotations'] = utils.order_vep_by_csq(variant['vep_annotations'])  # Adds major_consequence
            ordered_csqs = [x['major_consequence'] for x in variant['vep_annotations']]
            ordered_csqs = reduce(lambda x, y: ','.join([x, y]) if y not in x else x, ordered_csqs, '').split(',') # Close but not quite there
            consequences = defaultdict(lambda: defaultdict(list))
            for annotation in variant['vep_annotations']:
                annotation['HGVS'] = utils.get_proper_hgvs(annotation)
                consequences[annotation['major_consequence']][annotation['Gene']].append(annotation)
        assocset = str(db.list_association_sets()['name'][0])
        #vl = pandas.DataFrame([[chrom, pos, ref, alt]], columns = ['chrom','pos','ref','alt'])
        icdstats = db.get_association_data(association_set=assocset, chromosome =chrom, position=pos, pvalue_max = .01,  association_fields = ('pvalue', 'title', 'beta', 'odds_ratio', 'beta', 'se') )
        #icdstats = db.get_association_data(association_set=assocset, variant_list = vl, association_fields = ('pvalue', 'title', 'beta', 'odds_ratio', 'beta', 'se'))
        icdstats.query('ref == @refl and alt == @altl', inplace = True)
        icdstats = icdstats.to_dict('records')
        indexes = []
        seend = {}
        phef = db.get_phenotype_fields(association_set=str(db.list_association_sets()['name'][0]))
        for idx in range(0,len(icdstats)):
            # ICD10=T81/Name=Complications of procedures, not elsewhere classified/Chapter=T/L95OR=0.97/U95OR=2.04/OR=1.40/pvalue=0.0756463/l10pval=1.12/Case=1229
            item = icdstats[idx]
            icd10 = str(item['title'])
            item['Code'] = icd10
            if 'title' not in item or (numpy.isnan(item['odds_ratio']) and numpy.isnan(item['beta'])) or icd10 == "HC90":
                item['Name'] = 'NA'
                item['Group'] = 'NA'
                item['OR'] = 1
                item['LOR'] = 0
                item['L95OR'] = 1
                item['U95OR'] = 1
                item['pvalue'] = 1
                item['l10pval'] = 0
                item['Case'] = 'NA'
                item['SE'] = 0
                indexes.append(idx)
            else:
                sqz = str(phef[phef['title'] == icd10]['notes'].squeeze()).split(';')
                print(phef[phef['title'] == icd10])
                print(str(phef[phef['title'] == icd10]['notes'].squeeze()))
                item['Name'] = sqz[1].split('=')[1]
                item['Case'] = sqz[0].split('=')[1]
                item['Group'] = sqz[-1].split('=')[1]
               # item['Group'] = icd10 # default value
               # groups = ['RH', 'FH', 'HC', 'cancer', 'ADD', 'INI_FC', 'QT_FC', 'BIN_FC', 'INI', 'MED', 'BIN', 'BRMRI', 'BROADBIN', 'BROADQT']
               # for groupidx in range(0, len(groups)):
               #     if icd10.startswith(groups[groupidx]):
               #         item['Group'] = groups[groupidx]
               #         break
                if numpy.isnan(item['odds_ratio']):
                    item['OR'] = 1
                else:
                    item['OR'] = float(item['odds_ratio'])
                if numpy.isnan(item['beta']):
                    item['LOR'] = numpy.log(float(item['OR']))
                else:
                    item['LOR'] = float(item['beta'])
                item['L95OR'] = float(item['LOR']) - 1.96*float(item['se'])
                item['U95OR'] = float(item['LOR']) + 1.96*float(item['se'])
                item['pvalue'] = float(item['pvalue'])
                item['l10pval'] = -numpy.log10(float(item['pvalue']))
                item['SE'] = float(item['se'])
                if float(item['pvalue']) == 0:
#                    item['pvalue'] = numpy.finfo(float).eps
                    item['pvalue'] = 1e-250
                    item['pvalue'] = float(item['pvalue'])
                    item['l10pval'] = 250
                # item['Case'] = icd10info[0]['Case']
                item['log10pvalue'] = float(item['l10pval'])
                se =  float(item['se'])
                if item['l10pval'] < 1 or se >= .5 or (se >= .08 and item['OR'] == item['LOR']) or int(item['Case']) <= 100  or item['Code'] == "HC67" or icd10 in seend:
                    indexes.append(idx)
                seend[icd10] = icd10
        for index in sorted(indexes, reverse=True):
            del icdstats[index]
        print('Rendering variant: %s' % variant_str)
        print(variant)
        print(icdstats)
        print(consequences)
        print(namespace)
        return render_template(
            'variant.html',
            variant=variant,
            icdstats=icdstats,
            consequences=consequences,
            ordered_csqs=ordered_csqs,
            debug_message='',
            namespace=namespace,
            #plot_pval_data = [{'x': [1, 2, 3], 'y': [1,2, 3]}],
            pval_slider_max = variant_page_pval_slider_max(icdstats),
            plot_pval_data = variant_page_plot_pval_data(icdstats),
            plot_lor_data  = variant_page_plot_lor_data(icdstats),
        )
    except Exception as e:
        print('Failed on variant:', variant_str, '; Error=', traceback.format_exc())
        abort(404)

def variant_page_data_prep_sub(icdstats, sort_key='log10pvalue'):
    plot_d_raw = collections.defaultdict(list)
    keys = icdstats[0].keys()
    for key in keys:
        plot_d_raw[key] = np.array([x[key] for x in icdstats])
    plot_df = pandas.DataFrame(plot_d_raw).sort_values(
        by=['Group', sort_key], ascending=[True, False]
    )
    plot_d_dict = collections.defaultdict(collections.defaultdict)
    groups = sorted(set(plot_df['Group']))
    for group in groups:
        for key in keys:
            plot_d_dict[group][key] = list(plot_df[plot_df['Group'] == group][key])
    for group in groups:
        for key in ['OR', 'LOR', 'L95OR', 'U95OR', 'pvalue', 'SE', 'log10pvalue']:
            plot_d_dict[group][key] = [float(x) for x in plot_d_dict[group][key]]
    for group in groups:
        #error_bar = {'L95OR': -1, 'U95OR': 1}
        #for key in error_bar.keys():
        #    diff = np.array(plot_d_dict[group][key]) - np.array(plot_d_dict[group]['LOR'])
        #    plot_d_dict[group]['d{}'.format(key)] = [0 if np.isnan(x) else np.abs(x) for x in diff]
        plot_d_dict[group]['196SE'] = list( 1.96 * np.array(plot_d_dict[group]['SE']) )
    for group in groups:
        if group not in set(['Disease_outcome', 'cancer', 'Mental_health', 'Family_history','Others', 'Health_and_medical_history', 'Psychosocial_factors', 'Digestive_health']):
            beta_or_lor = 'BETA'
            beta_or_lor_val = plot_d_dict[group]['LOR']
            beta_or_lor_l95 = np.array(plot_d_dict[group]['LOR']) - 1.96 * np.array(plot_d_dict[group]['SE'])
            beta_or_lor_u95 = np.array(plot_d_dict[group]['LOR']) + 1.96 * np.array(plot_d_dict[group]['SE'])
        else:
            beta_or_lor = 'OR'
            beta_or_lor_val = np.exp(np.array(plot_d_dict[group]['LOR']))
            beta_or_lor_l95 = np.exp(np.array(plot_d_dict[group]['LOR']) - 1.96 * np.array(plot_d_dict[group]['SE']))
            beta_or_lor_u95 = np.exp(np.array(plot_d_dict[group]['LOR']) + 1.96 * np.array(plot_d_dict[group]['SE']))
        group_len = len(plot_d_dict[group]['Code'])
        plot_d_dict[group]['text'] = [
            '{}. Case: {}, P-value: {:.3e}, {} = {:.5f} (95% [{:.5f}, {:.5f}]), SE = {:.5f}'.format(
                ''.join([c if c != '_' else ' ' for c in x[0]]), x[1], x[2], x[3], x[4], x[5], x[6], x[7]
            ) for x in zip(
                plot_d_dict[group]['Name'],
                plot_d_dict[group]['Case'],
                plot_d_dict[group]['pvalue'],
                [beta_or_lor] * group_len,
                beta_or_lor_val,
                beta_or_lor_l95,
                beta_or_lor_u95,
                plot_d_dict[group]['SE'],
                #plot_d_dict[group]['L95OR'],
                #plot_d_dict[group]['U95OR'],
            )
        ]
    return plot_d_dict


def variant_page_pval_slider_max(icdstats):
#    return 100
    plot_d_dict = variant_page_data_prep_sub(icdstats, sort_key='log10pvalue')
    return np.max([np.max(v['log10pvalue']) for v in plot_d_dict.values()])


def variant_page_plot_pval_data(icdstats):
    plot_d_dict = variant_page_data_prep_sub(icdstats, sort_key='log10pvalue')
    groups = plot_d_dict.keys()
    plot_d = [{
        'x':    plot_d_dict[group]['Code'],
        'y':    plot_d_dict[group]['l10pval'],
        'text': plot_d_dict[group]['text'],
        'name': group,
        'type': 'scatter',
        'mode': 'markers',
        'marker': {'size': 16, },
        'hoverinfo':'x+text',
    } for group in groups]

    return plot_d


def variant_page_plot_lor_data(icdstats):
    plot_d_dict = variant_page_data_prep_sub(icdstats, sort_key='log10pvalue')
    groups = plot_d_dict.keys()
    plot_d = [{
        'x':    plot_d_dict[group]['Code'],
        'y':    plot_d_dict[group]['LOR'],
        'error_y': {
            'type'       : 'data',
            'symmetric'  : 'false',
            'array'      : plot_d_dict[group]['196SE'],
            'arrayminus' : plot_d_dict[group]['196SE'],
        },
        'text': plot_d_dict[group]['text'],
        'name': group,
        'type': 'scatter',
        'mode': 'markers',
        'marker': {'size': 16, },
        'hoverinfo':'x+text',
    } for group in groups]

    return plot_d


@app.route('/<namespace>/coding/<icd_str>')
def icd_page(namespace,icd_str):
#    if not check_credentials():
#        return redirect(url_for('login'))
    db = get_db(namespace)
        # icdlabel = str(icd_str)
        # .strip('ADD').strip('INI').strip('BRMRI')
        # Try several cutoffs to see if there are too many variants to
        # render.  Arbitrary 10k max.
        # passing = False
    cutoff = None
    icd = None
    for p in [.000001]:
        assocset = str(db.list_association_sets()['name'][0])
        df = db.get_phenotype_fields(association_set=assocset, include_pvalue_threshold=True)
        field_identifier = int(df[df['title'] == icd_str]['field_id'])
        pvalthr = p
       # pvalthr = float(df[df['title'] == icd_str]['pvalue_threshold'])
       # if np.isnan(pvalthr):
       #    pvalthr = p
        pvalthr = max(5e-7, pvalthr)
        icd = lookups.get_icd_variant_by_icd_id_pvalue(db, icd_str, field_identifier, pvalthr)
#        icd2 = icd
#        if icd2[~pandas.isnull(icd2.ld)].empty:
#            icd2 = icd2.query('log10pvalue > 7.5', inplace = False)
#        else:
#            icd2 = icd2.query('log10pvalue > 7.5 and ld', inplace = False)
   #     print(icd)
        icd = icd.to_dict(orient='records')
   #     print(icd)
#        icd2 = icd2.to_dict('records')
        if len(icd):
            cutoff = pvalthr
            break
    if icd is None or len(icd) == 0:
        icd = [{'Case': 'NA', 'Name': 'NA', 'icd': icd_str}]
    return render_template(
        'icd.html',
        icd=icd,
        icd2=icd,
        cutoff=cutoff,
        namespace=namespace)
#    except Exception as e:
#        print('Failed on icd:', icd_str, '; Error=', traceback.format_exc())
#        abort(404)



@app.route('/<namespace>/coding/dprs/<icd_str>')
def icddprs_page(namespace,icd_str):
    try:
        db = get_db(namespace)
        assocset = str(db.list_association_sets()['name'][0])
        df = db.get_phenotype_fields(association_set=assocset, include_pvalue_threshold=True)
        name = str(df[df['title'] == icd_str]['notes'].squeeze()).split(';')[1].split('=')[1]
        countn = str(df[df['title'] == icd_str]['notes'].squeeze()).split(';')[0].split('=')[1]
        cutoff = None
        icd = None
        print('Rendering ICD10: %s' % icd_str)
        if icd is None or len(icd) == 0:
            icd = [{'Case': countn, 'Name': name, 'icd': icd_str}]
        # print(icd_info)
        return render_template(
            'icddprs.html',
            icd=icd,
            namespace=namespace,
            cutoff=cutoff
            )
    except Exception as e:
        print('Failed on icd:', icd_str, '; Error=', traceback.format_exc())
        abort(404)



@app.route('/<namespace>/codinguncertain/<icd_str>')
def icduncertain_page(namespace,icd_str):
    db = get_db(namespace)
    try:
        cutoff = None
        icd = None

        for p in [.001]:
            icd = lookups.get_icd_variant_by_icd_id_pvalue_uncertain(db, icd_str, p)
            if len(icd):
                cutoff = p
                break
        print('Rendering ICD10: %s' % icd_str)

        if icd is None or len(icd) == 0:
            icd = [{'Case': 'NA', 'Name': 'NA', 'icd': icd_str}]

        return render_template(
            'icduncertain.html',
            icd=icd,
            namespace=namespace,
            cutoff=cutoff
            )
    except Exception as e:
        print('Failed on icd:', icd_str, '; Error=', traceback.format_exc())
        abort(404)



@app.route('/<namespace>/coding/gene/<icd_str>')
def codinggene_page(namespace,icd_str):
    db = get_db(namespace)
    try:
        cutoff = None
        icd = None
        for p in [.00001]:
            icd = lookups.get_icd_variant_by_icd_id_pvalue(db, icd_str, p)
            if len(icd):
                cutoff = p
                break
        print('Rendering ICD10: %s' % icd_str)

        if icd is None or len(icd) == 0:
            icd = [{'Case': 'NA', 'Name': 'NA', 'icd': icd_str}]
        # print(icd_info)

        return render_template(
            'icdgene.html',
            icd=icd,
            namespace=namespace,
            cutoff=cutoff
            )
    except Exception as e:
        print('Failed on icd:', icd_str, '; Error=', traceback.format_exc())
        abort(404)



@app.route('/<namespace>/coding/phenotype/<icd_str>')
def codingphenotype_page(namespace,icd_str):
    db = get_db(namespace)
    try:
        cutoff = None
        icd = None
        for p in [.00001]:
            assocset = str(db.list_association_sets()['name'][0])
            df = db.get_phenotype_fields(association_set=assocset, include_pvalue_threshold=True)
            field_identifier = int(df[df['title'] == icd_str]['field_id'])
            pvalthr = float(df[df['title'] == icd_str]['pvalue_threshold'])
            if np.isnan(pvalthr):
                pvalthr = p
            pvalthr = max(5e-7, pvalthr)
            icd = lookups.get_icd_variant_by_icd_id_pvalue(db, icd_str, field_identifier, pvalthr)
            if len(icd):
                cutoff = pvalthr
                break
        if icd is None or len(icd) == 0:
            icd = [{'Case': 'NA', 'Name': 'NA', 'icd': icd_str}]
        else:
            icd = icd.to_dict('records')
        return render_template(
            'phenotype.html',
            icd=icd,
            namespace=namespace,
            cutoff=cutoff
            )
    except Exception as e:
        print('Failed on icd:', icd_str, '; Error=', traceback.format_exc())
        abort(404)


@app.route('/<namespace>/target/1')
def target_page(namespace):
    #get_icd_variant_by_beta_pvalue(db, pvalue=0.00000005, betaabs = -.2)
    db = get_db(namespace)
    cutoff = None
    icd = None
    for p in [.000005]:
        assocset = str(db.list_association_sets()['name'][0])
        df = db.get_phenotype_fields(association_set=assocset, include_pvalue_threshold=True)
        pvalthr = p
        pvalthr = max(5e-7, pvalthr)
        icd = lookups.get_icd_variant_by_beta_pvalue(db, p, betaabs = -.1)
        icd.to_csv('tmp_target.csv',index=False)
        icd = icd.to_dict(orient='records')
        if len(icd):
            cutoff = pvalthr
            break
    if icd is None or len(icd) == 0:
        icd = [{'Case': 'NA', 'Name': 'NA', 'icd': icd_str}]
    return render_template(
        'target.html',
        icd = icd,
        namespace=namespace,
        cutoff=cutoff
    )


@app.route('/power')
def power_page():
    namespace = 'RIVAS_HG19'
    try:
        return render_template(
            'power.html',
            namespace = namespace
            )
    except Exception as e:
        print('Failed on protective scan Error=', traceback.format_exc())
        abort(404)

@app.route('/sex-effects')
def sex_page():
    namespace = 'RIVAS_HG19'
    try:
        return render_template(
            'sexeffects.html',
            namespace = namespace
            )
    except Exception as e:
        print('Failed on sex effects Error=', traceback.format_exc())
        abort(404)

@app.route('/mrpshiny')
def mrpshiny_page():
    namespace = 'RIVAS_HG19'
    try:
        return render_template(
            'mrpshiny.html',
            namespace = namespace
            )
    except Exception as e:
        print('Failed on MRP Shiny Error=', traceback.format_exc())
        abort(404)


@app.route('/snpnetcox')
def snpnetcox_page():
    namespace = 'RIVAS_HG19'
    try:
        return render_template(
            'snpnetcox.html',
            namespace = namespace
            )
    except Exception as e:
        print('Failed on MRP Shiny Error=', traceback.format_exc())
        abort(404)

@app.route('/cgauge')
def cguage_page():
    namespace = 'RIVAS_HG19'
    try:
        return render_template(
            'cgauge.html',
            namespace = namespace
            )
    except Exception as e:
        print('Failed on cgauge Error=', traceback.format_exc())
        abort(404)

@app.route('/biomarkers')
def biomarkers_page():
    namespace = 'RIVAS_HG19'
    try:
        return render_template(
            'biomarkers.html',
            namespace = namespace
            )
    except Exception as e:
        print('Failed on cgauge Error=', traceback.format_exc())
        abort(404)

@app.route('/<namespace>/snpnet/<icd_str>')
def snpnet_page(namespace, icd_str):
    namespace = 'RIVAS_HG19'
    db = get_db(namespace)
    try:
        cutoff = None
        icd = None
        for p in [.001]:
            assocset = str(db.list_association_sets()['name'][0])
            df = db.get_phenotype_fields(association_set=assocset, include_pvalue_threshold=True)
            field_identifier = int(df[df['title'] == icd_str]['field_id'])
            pvalthr = float(df[df['title'] == icd_str]['pvalue_threshold'])
            shortname = str(df[df['title'] == icd_str]['notes'].squeeze().split(';')[1].split('=')[1])
            casecnt = str(df[df['title'] == icd_str]['notes'].squeeze().split(';')[0].split('=')[1])
            pvalthr = max(5e-7, pvalthr)
            cuttoff = pvalthr
        icd = [{'Case': casecnt, 'Name': shortname, 'icd': icd_str}]

	return render_template(
            'snpnet.html',
            namespace=namespace,
            icd=icd,
            icd_str=icd_str,
            snpnet_plot='/static/PRS_map/{}.plot.png'.format(icd_str),
            snpnet_eval='/static/PRS_map/{}.eval.tsv'.format(icd_str)
        )
    except Exception as e:
        print('Failed on snpnet.html  Error=', traceback.format_exc())
        abort(404)

@app.route('/<namespace>/coding_breakdown/<icd_str>')
def coding_breakdown_page(namespace, icd_str):
    namespace = 'RIVAS_HG19'
#    if not check_credentials():
#        return redirect(url_for('login'))
    db = get_db(namespace)
    cutoff = None
    icd = None
    for p in [.000001]:
        assocset = str(db.list_association_sets()['name'][0])
        df = db.get_phenotype_fields(association_set=assocset, include_pvalue_threshold=True)
        field_identifier = int(df[df['title'] == icd_str]['field_id'])
        pvalthr = p
        pvalthr = max(5e-7, pvalthr)
        icd = lookups.get_icd_variant_by_icd_id_pvalue(db, icd_str, field_identifier, pvalthr)
        icd = icd.to_dict(orient='records')
        if len(icd):
            cutoff = pvalthr
            break
    if icd is None or len(icd) == 0:
        icd = [{'Case': 'NA', 'Name': 'NA', 'icd': icd_str}]
    try:
        return render_template(
            'upset.html',
            icd=icd,
            namespace=namespace
        )
    except Exception as e:
        print('Failed on upset.html  Error=', traceback.format_exc())
        abort(404)


@app.route('/degas')
def decomposition_dev_page():
        # this list file will be automatically updated by cron job that calls the following script:
        # https://github.com/rivas-lab/decomposition/blob/master/src/decomposition_dataset_list.sh
    with open('/biobankengine/app/static/decomposition/decomposition_datasets.lst') as f:
        dataset_list = f.read().splitlines()
    namespace = 'RIVAS_HG19'
    return render_template(
        'decomposition-dev.html',
        namespace=namespace,
        dataset_list = dataset_list,
    )


@app.route('/decomposition-internal/<dataset>')
def decomposition_internal_page(dataset):

    init_idx_pc  = 0
    init_idx_phe = 0
    init_idx_var = 0
    debug_str = 'debug'
    namespace = 'RIVAS_HG19'
    try:
        return render_template(
            'decomposition-internal.html',
            init_idx_pc  = init_idx_pc,
            init_idx_phe = init_idx_phe,
            init_idx_var = init_idx_var,
            dataset = dataset,
            namespace = namespace,
            debug_str = debug_str
        )

    except Exception as e:
        print('Unknown Error=', traceback.format_exc())
        abort(404)


@app.route('/decomposition-app')
def decomposition_app_page():
    abort(404)

def decomposition_app_page_null():
    db = get_db()

    try:
        with open('./static/decomposition/decomposition-app.lst') as f:
            dataset_list = f.read().splitlines()
        namespace = 'RIVAS_HG19'
        return render_template(
            'decomposition-app.html',
            dataset_list = dataset_list,
            namespace = namespace
        )

    except Exception as e:
        print('Unknown Error=', traceback.format_exc())
        abort(404)



@app.route('/decomposition/<dataset>')
def decomposition_page(dataset):

    if(dataset == 'PTVs'):
        init_idx_pc  = 6
        init_idx_phe = 429
        init_idx_var = 8503
    else:
        init_idx_pc  = 0
        init_idx_phe = 0
        init_idx_var = 0

    debug_str = 'debug'
    namespace = 'RIVAS_HG19'
    try:
        if(dataset == "20170930_EMBL-Stanford_coding-nonMHC_z"):
            return render_template(
                'decomposition-20170930.html',
                init_idx_pc  = init_idx_pc,
                init_idx_phe = init_idx_phe,
                init_idx_var = init_idx_var,
                dataset = dataset,
                namespace = namespace,
                debug_str = debug_str
            )
        elif(dataset[:8] == "20171011"):
            return render_template(
                'decomposition-20171011.html',
                init_idx_pc  = init_idx_pc,
                init_idx_phe = init_idx_phe,
                init_idx_var = init_idx_var,
                dataset = dataset,
                namespace = namespace,
                debug_str = debug_str
                )
        else:
            return render_template(
                'decomposition.html',
                namesspace = namespace,
                init_idx_pc  = init_idx_pc,
                init_idx_phe = init_idx_phe,
                init_idx_var = init_idx_var,
                dataset = dataset,
                debug_str = debug_str
            )

    except Exception as e:
        print('Unknown Error=', traceback.format_exc())
        abort(404)


@app.route('/hla-assoc')
def hla_assoc_page():
    namespace = 'RIVAS_HG19'
    try:
        return render_template('hla-assoc.html',
                               namespace = namespace)

    except Exception as e:
        print('Unknown Error=', traceback.format_exc())
        abort(404)


@app.route('/prs')
def prs_page():
    # Please look at the documentation at PRS the repository
    # https://github.com/rivas-lab/PRS/tree/master/notebook/20210115_GBE_data_prep/GBE_page_src
    try:
        namespace = 'RIVAS_HG19'
        if request.method == 'POST':
            namespace = request.form['functionassocset']

        # read trait  list table
        trait_list_f='/biobankengine/app/static/PRS_map/traits.tsv'
        table_cols=['Trait group', 'Trait', 'Family', 'Geno', 'Covars', 'Full', 'delta', '# variants', 'p (WB)', 'significant?']
        table_cols_select=['Trait group', 'Family', 'significant?']

        df = pandas.read_csv(trait_list_f, sep='\t')
        df['trait'] = ['<a href="/RIVAS_HG19/snpnet/{}">{}</a>'.format(x[0], x[1]) for x in zip(df['trait'], df['trait_name'])]
        df = df.drop('trait_name',axis=1)
        for col in ['#trait_category']:
                # format string
            df[col] = df[col].map(lambda x: str(x).replace('_', ' '))
        for col in ['WB_test_P']:
                # format to scientific notation
            df[col] = df[col].map(lambda x: '{:0.2e}'.format(x))
        for col in ['geno', 'covar', 'geno_covar', 'geno_delta']:
                # format digits
            df[col] = df[col].map(lambda x: str(round(x, 2)))

        # generate HTML string
        table_prs_trait_list_tbody_str=''.join(['<tr>{}</tr>'.format(
                ''.join(['<td>{}</td>'.format(x) for x in df.iloc[row]])
        ) for row in range(df.shape[0])])

        return render_template(
            'prs.html',
            namespace = namespace,
            table_prs_trait_list_cols        = table_cols,
            table_prs_trait_list_col_len     = len(table_cols),
            table_prs_trait_list_cols_select = table_cols_select,
            table_prs_trait_list_tbody_str   = table_prs_trait_list_tbody_str
	    )

    except Exception as e:
        print('Unknown Error=', traceback.format_exc())
        abort(404)


@app.route('/<namespace>/mrpgene/<freq>')
def mrp_page(namespace,freq):
    try:
	# read trait  list table
        if freq == "rare":
	    trait_list_f='/biobankengine/app/static/mrpgene/' + namespace + '/' + namespace + '_ultrarare_gene_results.tsv'
        else:
            trait_list_f='/biobankengine/app/static/mrpgene/' + namespace + '/' + namespace + '_gene_results.tsv'
	table_cols=['GBE_ID','GBE_short_name','pops','num_pops','gene','n_pav','n_ptv','l10BF_IEM_pav','l10BF_SEM_pav','l10BF_SEM_ptv','GBE_link','MRPMM_pav','MRPMM_ptv']
        df = pandas.read_csv(trait_list_f, sep='\t')
        df['gene2'] = df['gene']
        df['gene'] = ['<a href="/' + namespace + '/gene/{}">{}</a>'.format(x[0], x[1]) for x in zip(df['gene2'], df['gene2'])]
	df['trait'] = ['<a href="/' + namespace + '/coding/{}">{}</a>'.format(x[0], x[1]) for x in zip(df['GBE_ID'], df['GBE_short_name'])]
        df['mrpmm_pav'] = ['<a href="/' + namespace + '/mrpmm_pav/{}/{}">{}</a>'.format(x[0],x[1], x[2]) for x in zip(df['gene2'], df['GBE_ID'], df['GBE_short_name'])]
        df['mrpmm_ptv'] = ['<a href="/' + namespace + '/mrpmm_ptv/{}/{}">{}</a>'.format(x[0],x[1], x[2]) for x in zip(df['gene2'], df['GBE_ID'], df['GBE_short_name'])]
        del df['gene2']
	# generate HTML string
	table_mrp_trait_list_tbody_str=''.join(['<tr>{}</tr>'.format(
    	    ''.join(['<td>{}</td>'.format(x) for x in df.iloc[row]])
	) for row in range(df.shape[0])])
        if namespace == "RIVAS_HG38":
            return render_template(
                'mrp200k.html',
		namespace = namespace,
        table_mrp_trait_list_cols        = table_cols,
		table_mrp_trait_list_col_len     = len(table_cols),
#		table_mrp_trait_list_cols_select = table_cols_select,
		table_mrp_trait_list_tbody_str   = table_mrp_trait_list_tbody_str
	    )
        else:
            return render_template(
                'mrparray.html',
		namespace = namespace,
        table_mrp_trait_list_cols        = table_cols,
		table_mrp_trait_list_col_len     = len(table_cols),
#		table_mrp_trait_list_cols_select = table_cols_select,
		table_mrp_trait_list_tbody_str   = table_mrp_trait_list_tbody_str
	    )


    except Exception as e:
        print('Unknown Error=', traceback.format_exc())
        abort(404)

@app.route('/dprs')
def dprs_page():
    try:
        namespace = 'RIVAS_HG19'
        if request.method == 'POST':
            namespace = request.form['functionassocset']
        return render_template('dprs.html',
                               namespace = namespace)

    except Exception as e:
        print('Unknown Error=', traceback.format_exc())
        abort(404)


def get_icd_variant_table(icd, namespace, p):
    db = get_db(namespace)
    significant_variant_ids = {}
    significant_variants = lookups.get_icd_significant(db, icd, p)
    for v in significant_variants:
        significant_variant_ids[v['xpos']] = {}
        significant_variant_ids[v['xpos']]['pvalue'] = v['pvalue']
        significant_variant_ids[v['xpos']]['or'] = v['or_val']
    variants = lookups.get_variants_by_id(db, significant_variant_ids.keys())
    for v in variants:
        genes = []
        symbols = []
        for a in v['vep_annotations']:
            if not a['Gene'] in genes:
                genes.append(a['Gene'])
                symbols.append(a['SYMBOL'])


        v['gene_name'] = ",".join(genes[0:3])
        v['gene_symbol'] = ",".join(symbols[0:3])
        v['this_icd'] = icd
        v['this_pvalue'] = significant_variant_ids[v['xpos']]['pvalue']
        v['this_or'] = significant_variant_ids[v['xpos']]['or']
    return variants


def get_prot_variant_table(lor, namespace, p):
    db = get_db(namespace)
    significant_variant_ids = {}
    significant_variants = lookups.get_significant_prot(db,p, lor)
    for v in significant_variants:
        significant_variant_ids[v['xpos']] = {}
        significant_variant_ids[v['xpos']]['pvalue'] = v['stats'][0]['pvalue']
        significant_variant_ids[v['xpos']]['or'] = v['stats'][0]['or']
    variants = lookups.get_variants_by_id(db, significant_variant_ids.keys())
    for v in variants:
        genes = []
        symbols = []
        for a in v['vep_annotations']:
            if not a['Gene'] in genes:
                genes.append(a['Gene'])
                symbols.append(a['SYMBOL'])


        v['gene_name'] = ",".join(genes[0:3])
        v['gene_symbol'] = ",".join(symbols[0:3])
        v['this_pvalue'] = significant_variant_ids[v['xpos']]['pvalue']
        v['this_or'] = significant_variant_ids[v['xpos']]['or']
    return variants



@app.route('/<namespace>/intensity/<affy_str>')
def intensity_page(namespace, affy_str):
    db = get_db(namespace)
    try:
        n_UKBL = 11
        n_UKBB = 95


        print('Rendering Intensity page: %s' % affy_str)
        #variant = lookups.get_icd_affyid(db, affy_str)
        # print(variant)
        return render_template(
            'intensity.html',
            affy=affy_str,
         #   variant=variant,
            namespace = namespace,
            UKBL_idx = [1 + x for x in range(n_UKBL)],
            UKBB_idx = [1 + x for x in range(n_UKBB)]
            )
    except Exception as e:
        print('Failed on affy id:', affy_str, '; Error=', traceback.format_exc())
        abort(404)

@app.route('/<namespace>/gene/<gene_id>', methods=['GET', 'POST'])
def gene_page(namespace, gene_id):
    if gene_id in []:
        return open(os.path.join(GENE_CACHE_DIR, '{}.html'.format(gene_id))).read()
    else:
        if request.method == 'POST':
            if 'submit_phenos' in request.form.keys():
                functionphen = request.form.getlist('phenotypes[]')
                phenidarr = []
                if len(functionphen) == 0:
                    phenidarr.append("ALL")
                for phenname in functionphen:
                    phenidarr.append(str(phenname))
            else:
                phenidarr.append('ALL')
            if len(phenidarr) == 0:
                abort(404)
            elif 'ALL' in phenidarr:
                return get_gene_page_content_all(namespace, gene_id, phenidarr)
            else:
                return get_gene_page_content_all(namespace, gene_id, phenidarr)
        else:
            phenidarr = []
            phenidarr.append("ALL")
            return get_gene_page_content_all(namespace, gene_id, phenidarr)



@app.route('/<namespace>/mrpmm_pav/<gene_id>/<icd_str>')
def mrpmm_pav(namespace, gene_id, icd_str):
    l10bfarr = []
    with open(glob.glob("/biobankengine/app/static/mrpmm/" + namespace + "_mrpmm/" + gene_id + "_" + icd_str + "*pav*.mcmc.bic.aic")[0],"r") as inFile:
        finr = inFile.readlines()
       # l10bf = {}
        for line in finr[1:]:
            line = line.rstrip()
            line = line.split()
            #            l10bf[line[0]] = "{0:.3g}".format(float(line[1]))
            l10bf = "{0:.3g}".format(float(line[len(line) - 1]))
            l10bfarr.append(float(l10bf))
    print(l10bfarr)
    idxmax = 1
    l10bfarrtmp = 0
    for i in range(0, len(l10bfarr)):
        if l10bfarr[i] >= 2.5 and l10bfarr[i] > l10bfarrtmp + 2:
            idxmax = i + 1
            l10bfarrtmp = l10bfarr[i]
        else:
            pass
#    idxmax = l10bfarr.index(max(l10bfarr)) + 1
    #l10bf = max(l10bfarr)
    l10bf = l10bfarrtmp
    with open(glob.glob("/biobankengine/app/static/mrpmm/" + namespace + "_mrpmm/"+ gene_id + "_" + icd_str + "_*pav*" + str(idxmax) + ".mcmc.posteriors")[0], "r") as inFile:
        probabilities = []
        admixture_data = []
        variants = []
        nullmembership = []
        finr = inFile.readlines()
        for line in finr[1:]:
            line = line.rstrip()
            var_info = line.split("\t")
            varid = var_info[0]
            var = var_info[4]
            variant = {}
            varids = varid.split(':')
            varidpage = varids[0] + ":" + varids[1] + "-" + varids[2] + "-" + varids[3]
            variant["variant"] = varidpage
            variant["varid"] = var
            for j in range(5, len(var_info)):
                variant["%s" % (j-4)] = var_info[j]
            if float(var_info[5]) >= .2:
                        continue
            nullmembership.append(float(var_info[5]))
            variants.append(variant)
        idxnewarr = [b[0] for b in sorted(enumerate(nullmembership),key=lambda i:i[1], reverse = True)]
        variants = [variants[idxnewarr[i]] for i in range(0,len(idxnewarr))]
        admixture_data.append(variants)
    with open(glob.glob("/biobankengine/app/static/mrpmm/" + namespace +"_mrpmm/"+ gene_id + "_"  + icd_str + "_*pav*" + str(idxmax) + ".mcmc.gene.posteriors")[0], "r") as inFile:
        admixture_datagene = []
        genes = []
        finr = inFile.readlines()
        for line in finr[1:]:
            line = line.rstrip()
            gene_info = line.split("\t")
            geneid = gene_info[0]
            gene = {}
            gene["gene"] = geneid
            for j in range(1, int((len(gene_info)-1)/3)+1):
                gene["%s" % (j-1)] = gene_info[j]
            genes.append(gene)
        admixture_datagene.append(genes)
    scalefactor = 1
    scales = []
    with open(glob.glob("/biobankengine/app/static/mrpmm/" + namespace + "_mrpmm/" + gene_id + "_" + icd_str + "_*pav*" + str(idxmax) + ".mcmc.scale")[0], 'r') as inFile:
        finr = inFile.readlines()
        for line in finr[1:]:
            line = line.rstrip()
            line = line.split()
            scaled = float(line[2])
            scales.append(scaled)
    scalefactor = numpy.median(scales)
    with open(glob.glob("/biobankengine/app/static/mrpmm/" + namespace + "_mrpmm/"+ gene_id + "_" + icd_str + "_*pav*" + str(idxmax) + ".mcmc.bc")[0], "r") as inFile:
        cluster_data = []
        headerarr = []
        inFiler = inFile.readlines()
        header = inFiler[0]
        header = header.split('\t')
        for j in range(1,len(header)):
            if j % 3 == 1:
                headeritem = header[j].decode('utf8').replace(" ","")
                headeritem = headeritem.rstrip(header[j][-3:])
                headerarr.append(headeritem)
        for line in inFiler[1:]:
            line = line.rstrip()
            cluster = []
            clusterl95 = []
            clusteru95 = []
            clust_info = line.split("\t")
            cluster_num = clust_info[0]
            for j in range(1, len(clust_info)):
                if j % 3 == 1:
                    tmp = []
                    headeritem = header[j].decode('utf8').replace(" ","")
                    headeritem = headeritem.rstrip(header[j][-3:])
                    tmp.append(headeritem.encode("ascii"))
                    tmp.append(min(float(clust_info[j]),scalefactor*float(clust_info[j])))
                if j % 3 == 2:
                    tmp.append(min(float(clust_info[j]),scalefactor*float(clust_info[j])))
                if j % 3 == 0:
                    tmp.append(min(float(clust_info[j]),scalefactor*float(clust_info[j])))
                    cluster.append(tmp)
            cluster_data.append([cluster_num, cluster])
    fdrf = open(glob.glob("/biobankengine/app/static/mrpmm/" + namespace + "_mrpmm/"+ gene_id + "_" + icd_str +"_*pav*" + str(idxmax) + ".fdr")[0], "r").readlines()
    fdrnum = fdrf[0].rstrip().split()[0]
    fdr_data2 = []
    for line in fdrf[1:]:
        line = line.rstrip()
        line = line.split()
        fdr_data2.append(line[0])
    try:
        t = render_template(
            'mrp.html',
            namespace=namespace,
            plot_data=admixture_data,
            gene_data=admixture_datagene,
            num_figs=cluster_data,
            log10_bf=l10bf,
            cluster_num=str(int(cluster_num)+1),
            fdr=fdrnum,
            fdr_data=fdr_data2
   #         phenid_arr=headerarr
            )
        return t
    except Exception as e:
        print('Failed: %s' % e)
        abort(404)


@app.route('/<namespace>/mrpmm_ptv/<gene_id>/<icd_str>')
def mrpmm_ptv(namespace, gene_id, icd_str):
    l10bfarr = []
    with open(glob.glob("/biobankengine/app/static/mrpmm/" + namespace +"_mrpmm/" + gene_id + "_" + icd_str + "*ptv*.mcmc.bic.aic")[0],"r") as inFile:
        finr = inFile.readlines()
       # l10bf = {}
        for line in finr[1:]:
            line = line.rstrip()
            line = line.split()
            #            l10bf[line[0]] = "{0:.3g}".format(float(line[1]))
            l10bf = "{0:.3g}".format(float(line[len(line) - 1]))
            l10bfarr.append(float(l10bf))
    print(l10bfarr)
    idxmax = l10bfarr.index(max(l10bfarr)) + 1
    l10bf = max(l10bfarr)
    with open(glob.glob("/biobankengine/app/static/mrpmm/" + namespace +"_mrpmm/"+ gene_id + "_" + icd_str + "_*ptv*" + str(idxmax) + ".mcmc.posteriors")[0], "r") as inFile:
        probabilities = []
        admixture_data = []
        variants = []
        nullmembership = []
        finr = inFile.readlines()
        for line in finr[1:]:
            line = line.rstrip()
            var_info = line.split("\t")
            varid = var_info[0]
            var = var_info[4]
            variant = {}
            varids = varid.split(':')
            varidpage = varids[0] + ":" + varids[1] + "-" + varids[2] + "-" + varids[3]
            variant["variant"] = varidpage
            variant["varid"] = var
            for j in range(5, len(var_info)):
                variant["%s" % (j-4)] = var_info[j]
            if float(var_info[5]) >= .2:
                        continue
            nullmembership.append(float(var_info[5]))
            variants.append(variant)
        idxnewarr = [b[0] for b in sorted(enumerate(nullmembership),key=lambda i:i[1], reverse = True)]
        variants = [variants[idxnewarr[i]] for i in range(0,len(idxnewarr))]
        admixture_data.append(variants)
    with open(glob.glob("/biobankengine/app/static/mrpmm/" + namespace +"_mrpmm/"+ gene_id + "_"  + icd_str + "_*ptv*" + str(idxmax) + ".mcmc.gene.posteriors")[0], "r") as inFile:
        admixture_datagene = []
        genes = []
        finr = inFile.readlines()
        for line in finr[1:]:
            line = line.rstrip()
            gene_info = line.split("\t")
            geneid = gene_info[0]
            gene = {}
            gene["gene"] = geneid
            for j in range(1, int((len(gene_info)-1)/3)+1):
                gene["%s" % (j-1)] = gene_info[j]
            genes.append(gene)
        admixture_datagene.append(genes)
    scalefactor = 1
    scales = []
    with open(glob.glob("/biobankengine/app/static/mrpmm/" + namespace +"_mrpmm/" + gene_id + "_" + icd_str + "_*ptv*" + str(idxmax) + ".mcmc.scale")[0], 'r') as inFile:
        finr = inFile.readlines()
        for line in finr[1:]:
            line = line.rstrip()
            line = line.split()
            scaled = float(line[2])
            scales.append(scaled)
    scalefactor = numpy.median(scales)
    with open(glob.glob("/biobankengine/app/static/mrpmm/" + namespace +"_mrpmm/"+ gene_id + "_" + icd_str + "_*ptv*" + str(idxmax) + ".mcmc.bc")[0], "r") as inFile:
        cluster_data = []
        headerarr = []
        inFiler = inFile.readlines()
        header = inFiler[0]
        header = header.split('\t')
        for j in range(1,len(header)):
            if j % 3 == 1:
                headeritem = header[j].decode('utf8').replace(" ","")
                headeritem = headeritem.rstrip(header[j][-3:])
                headerarr.append(headeritem)
        for line in inFiler[1:]:
            line = line.rstrip()
            cluster = []
            clusterl95 = []
            clusteru95 = []
            clust_info = line.split("\t")
            cluster_num = clust_info[0]
            for j in range(1, len(clust_info)):
                if j % 3 == 1:
                    tmp = []
                    headeritem = header[j].decode('utf8').replace(" ","")
                    headeritem = headeritem.rstrip(header[j][-3:])
                    tmp.append(headeritem.encode("ascii"))
                    tmp.append(min(float(clust_info[j]),scalefactor*float(clust_info[j])))
                if j % 3 == 2:
                    tmp.append(min(float(clust_info[j]),scalefactor*float(clust_info[j])))
                if j % 3 == 0:
                    tmp.append(min(float(clust_info[j]),scalefactor*float(clust_info[j])))
                    cluster.append(tmp)
            cluster_data.append([cluster_num, cluster])
    fdrf = open(glob.glob("/biobankengine/app/static/mrpmm/" + namespace +"_mrpmm/"+ gene_id + "_" + icd_str +"_*ptv*" + str(idxmax) + ".fdr")[0], "r").readlines()
    fdrnum = fdrf[0].rstrip().split()[0]
    fdr_data2 = []
    for line in fdrf[1:]:
        line = line.rstrip()
        line = line.split()
        fdr_data2.append(line[0])
    try:
        t = render_template(
            'mrp.html',
            namespace=namespace,
            plot_data=admixture_data,
            gene_data=admixture_datagene,
            num_figs=cluster_data,
            log10_bf=l10bf,
            cluster_num=str(int(cluster_num)+1),
            fdr=fdrnum,
            fdr_data=fdr_data2
   #         phenid_arr=headerarr
            )
        return t
    except Exception as e:
        print('Failed: %s' % e)
        abort(404)


def get_gene_page_content_all(namespace, gene_id, phens):
    db = get_db(namespace)
    gene_id = str(gene_id)
    if True:
        t = None
        if t is None:
            geneidarr = []
            geneidarr.append(gene_id)
            genedf = db.get_genes(gene_name=geneidarr, exact_match = True)
            if genedf.empty:
                abort(404)
            gene_idx = genedf['gene_id']
            assocset = str(db.list_association_sets()['name'][0])
            phenidarr = []
            if 'ALL' in phens:
                dfphen = db.get_phenotype_fields(association_set=assocset)
            else:
                dfphen = db.get_phenotype_fields(association_set=assocset)
                for phetitle in phens:
                    field_identifier = int(dfphen[dfphen['title'] == phetitle]['field_id'])
                    phenidarr.append(field_identifier)
            dfv = db.get_variants(gene_name=gene_id, pad=100, variant_fields = ['ref', 'alt', 'rsid', 'gene', 'consequence', 'impact', 'hgvsp', 'lof', 'miss', 'miss_bileve', 'miss_wcsg', 'hwep', 'maf', 'ld', 'wcsg_only', 'bileve_only', 'missingness', 'hwe', 'mcpi', 'gnomad_filter', 'mgi', 'mgi_notes', 'all_filters', 'annotations', 'variant_identity', 'major_consequence', 'category', 'gene_name', 'gene_symbol', 'HGVSp', 'HGVSc', 'ukbb_freq'])
           # coverage_stats = lookups.get_coverage_for_transcript(
           #     db,
           #     transcript['xstart'] - EXON_PADDING,
           #     transcript['xstop'] + EXON_PADDING)
            if 'ALL' in phens:
                df = db.get_association_data(association_set=assocset, gene_name=gene_id, pad=500, pvalue_max=.00005)
                sefilter = .2
                columns = ['pvalue','title','beta','all_filters','se']
                ndf = dfv.apply(lambda row: pandas.Series([1,"NA",0,0,"NA"], index = columns) if df[df['variant_identity'] == row['variant_identity']]['pvalue'].empty else (pandas.Series([numpy.log(x) if i == 2 else x for i, x in enumerate(df.iloc[df[df['variant_identity'] == row['variant_identity']]['pvalue'].idxmin(),:][['pvalue','title','odds_ratio', 'all_filters','se']])], index = columns) if numpy.isnan(df.iloc[df[df['variant_identity'] == row['variant_identity']]['pvalue'].idxmin(),:]['beta']) and df[df['variant_identity'] == row['variant_identity']].query('se <= .2')['pvalue'].empty else (df.iloc[df[df['variant_identity'] == row['variant_identity']]['pvalue'].idxmin(),:][['pvalue','title','beta', 'all_filters','se']] if df[df['variant_identity'] == row['variant_identity']].query('se <= .2')['pvalue'].empty else (pandas.Series([numpy.log(x) if i == 2 else x for i, x in enumerate(df.iloc[df[df['variant_identity'] == row['variant_identity']].query('se <= .2')['pvalue'].idxmin(),:][['pvalue','title','odds_ratio', 'all_filters','se']])], index = columns) if numpy.isnan(df.iloc[df[df['variant_identity'] == row['variant_identity']].query('se <= .2')['pvalue'].idxmin(),:]['beta']) else ( df.iloc[df[df['variant_identity'] == row['variant_identity']].query('se <= .2')['pvalue'].idxmin(),:][['pvalue','title','beta', 'all_filters','se']])))), axis = 1)
                ndf.columns = ['minicd_p','minicd_icd','minicd_lor', 'filter', 'se']
                dfv = pandas.concat([dfv, ndf], axis = 1 , sort = False)
            else:
                df = db.get_association_data(association_set=assocset, gene_name=gene_id, pad=500, pvalue_max=1, field_id = phenidarr)
                sefilter = .2
                columns = ['pvalue','title','beta','all_filters','se']
                ndf = dfv.apply(lambda row: pandas.Series([1,"NA",0,0,"NA"], index = columns) if df[df['variant_identity'] == row['variant_identity']]['pvalue'].empty else (pandas.Series([numpy.log(x) if i == 2 else x for i, x in enumerate(df.iloc[df[df['variant_identity'] == row['variant_identity']]['pvalue'].idxmin(),:][['pvalue','title','odds_ratio', 'all_filters','se']])], index = columns) if numpy.isnan(df.iloc[df[df['variant_identity'] == row['variant_identity']]['pvalue'].idxmin(),:]['beta']) and df[df['variant_identity'] == row['variant_identity']].query('se <= .2')['pvalue'].empty else (df.iloc[df[df['variant_identity'] == row['variant_identity']]['pvalue'].idxmin(),:][['pvalue','title','beta', 'all_filters','se']] if df[df['variant_identity'] == row['variant_identity']].query('se <= .2')['pvalue'].empty else (pandas.Series([numpy.log(x) if i == 2 else x for i, x in enumerate(df.iloc[df[df['variant_identity'] == row['variant_identity']].query('se <= .2')['pvalue'].idxmin(),:][['pvalue','title','odds_ratio', 'all_filters','se']])], index = columns) if numpy.isnan(df.iloc[df[df['variant_identity'] == row['variant_identity']].query('se <= .2')['pvalue'].idxmin(),:]['beta']) else ( df.iloc[df[df['variant_identity'] == row['variant_identity']].query('se <= .2')['pvalue'].idxmin(),:][['pvalue','title','beta', 'all_filters','se']])))), axis = 1)
                ndf.columns = ['minicd_p','minicd_icd','minicd_lor', 'filter', 'se']
                dfv = pandas.concat([dfv, ndf], axis = 1 , sort = False)
            dfv['shortname'] = dfv.apply(lambda row: 'NA' if (row['minicd_icd'] == "NA" or len(str(dfphen[dfphen['title'] == row['minicd_icd']]['notes'].squeeze()).split(';')) <=  1) else str(dfphen[dfphen['title'] == row['minicd_icd']]['notes'].squeeze()).split(';')[1].split('=')[1], axis = 1)
            dfv['shortname'] = dfv.apply(lambda row: '_'.join(row['shortname'].split()), axis = 1)
            dfv['minpval'] = dfv['minicd_p']
            dfv['minor'] = dfv['minicd_lor']
            dfv['odds_ratio'] = dfv['minor']
            dfv['minpval'] = dfv.apply(lambda row: 1e-300 if row['minicd_p'] == 0 else row['minicd_p'], axis = 1)
            dfv['log10pvalue'] = dfv.apply(lambda row: -numpy.log10(row['minpval']), axis = 1)
            dfv['filter'] = dfv.apply(lambda row: 'PASS' if (row['filter'] == 0 and row['se'] != "NA" and row['se'] <= sefilter) else ("PASS" if row['se'] <= sefilter and row['se'] != "NA" else ('SE' if row['se'] > sefilter and row['se'] != "NA" else 'FAIL')), axis = 1)
            dfv['ukbb_freq'] = dfv['maf']
            transcripts_in_gene = []
            transcript_id = str(genedf['canonical_transcript_eid'][0])
         #   transcript = lookups.add_xpos(
         #       lookups.cast_pos_info(
         #           dict(zip(lookups.TRANSCRIPT_INFO_KEYS,
         #                    gene['c_transcript_info'].split(':')))))
         #   transcript['transcript_id'] = transcript_id
         #   transcript['exons'] = [
         #       lookups.cast_pos_info(
         #           dict(zip(lookups.EXON_INFO_KEYS, s.split(':'))))
         #       for s in gene['exon_info'].split(';')
         #       if s]
#            variants_in_transcript = lookups.get_variants_by_transcript_idx(
#                db, transcript_idx, transcript_id)
            #coverage_stats = dfv[['chrom','pos','filter','category','odds_ratio', 'log10pvalue']]
            #coverage_stats = coverage_stats.to_dict('records')
            #add_transcript_coordinate_to_variants(
            #    variants_in_transcript, transcript)
            variants_in_gene = dfv.to_dict('records')
            coverage_stats = variants_in_gene
            variants_in_transcript = variants_in_gene
            genedf = genedf.to_dict('records')[0]
            gbuild = namespace.split('_')[1].lower()
            t = render_template(
                'gene.html',
                gene=genedf,
                variants_in_gene=variants_in_gene,
                coverage_stats=coverage_stats,
                variants_in_transcript=variants_in_transcript,
                 gbuild=gbuild,
                namespace=namespace)
        print('Rendering gene: %s' % gene_idx)
        return t





#### Compare Betas for two phenotypes

@app.route('/<namespace>/comparebeta', methods=['GET', 'POST'])
def compare_page(namespace):
    if "pass" in []:
        return open(os.path.join(GENE_CACHE_DIR, '{}.html'.format(gene_id))).read()
    else:
        if request.method == 'POST':
            if 'submit_phenos' in request.form.keys():
                functionphen = request.form.getlist('phenotypes[]')
                phenidarr = []
                if len(functionphen) == 0:
                    abort(404)
                for phenname in functionphen:
                    phenidarr.append(str(phenname))
                if len(phenidarr) != 2:
                    abort(404)
            else:
                phenidarr = ['INI10030870','INI10030780']
            if len(phenidarr) == 0:
                abort(404)
            return get_compare_page_content_all(namespace, phenidarr)
        else:
            phenidarr = ['INI10030870','INI10030780']
            return get_compare_page_content_all(namespace, phenidarr)
#            abort(404)

def get_compare_page_content_all(namespace, phens):
    db = get_db(namespace)
    if True:
        t = None
        if t is None:
            assocset = str(db.list_association_sets()['name'][0])
            phenidarr = []
            if 'ALL' in phens:
                dfphen = db.get_phenotype_fields(association_set=assocset)
            else:
                dfphen = db.get_phenotype_fields(association_set=assocset)
                for phetitle in phens:
                    field_identifier = int(dfphen[dfphen['title'] == phetitle]['field_id'])
                    phenidarr.append(field_identifier)
            if 'PASS' in phens:
                a = 1
            else:
                df1 = db.get_association_data(association_set=assocset, pvalue_max=.00000005, field_id = phenidarr[0], variant_fields = ['ref','alt'])
                df2 = db.get_association_data(association_set=assocset, pvalue_max=.00000005, field_id = phenidarr[1], variant_fields = ['ref','alt'])
                res1 = df1[['chrom', 'pos', 'ref', 'alt']]
                res2 = df2[['chrom', 'pos', 'ref', 'alt']]
                resf = pandas.concat([res1, res2], axis = 0 , sort = False)
                df1 = db.get_association_data(association_set=assocset, field_id = phenidarr[0], variant_list = resf, variant_fields = ['chrom', 'pos', 'ref', 'alt', 'rsid', 'gene', 'consequence',  'ld',  'all_filters', 'variant_identity', 'major_consequence', 'category', 'gene_symbol', 'HGVSp', 'HGVSc', 'ukbb_freq'])
                df2 = db.get_association_data(association_set=assocset, field_id = phenidarr[1], variant_list = resf,  variant_fields = ['chrom', 'pos', 'ref', 'alt', 'rsid', 'gene', 'consequence',  'ld',  'all_filters', 'variant_identity', 'major_consequence', 'category', 'gene_symbol', 'HGVSp', 'HGVSc', 'ukbb_freq'])
                result = pandas.merge(df1, df2, left_on='variant_identity', right_on='variant_identity', how='outer')
                sefilter = .2
                columns = ['pvalue1','pvalue2', 'title1','title2', 'beta1','beta2', 'all_filters','se1', 'se2']
                ndf = result.apply(lambda row: pandas.Series([row['pvalue_x'] if not pandas.isnull(row['pvalue_x']) else 1, row['pvalue_y'] if not pandas.isnull(row['pvalue_y']) else 1, row['title_x'] if not pandas.isnull(row['pvalue_x']) else "NA", row['title_y'] if not pandas.isnull(row['title_y']) else "NA", numpy.log(row['odds_ratio_x']) if pandas.isnull(row['beta_x']) and not pandas.isnull(row['odds_ratio_x']) else (row['beta_x'] if not pandas.isnull(row['beta_x']) else 0), numpy.log(row['odds_ratio_y']) if pandas.isnull(row['beta_y']) and not pandas.isnull(row['odds_ratio_y']) else (row['beta_y'] if not pandas.isnull(row['beta_y']) else 0), 'MHC' if (int(row['chrom_x']) == 6 and int(row['pos_x']) >= 27000000 and int(row['pos_x']) <= 33000000) else ('PASS' if (row['se_x'] <= .2  or pandas.isnull(row['se_x'])) and (row['se_y'] <= .2 or pandas.isnull(row['se_y'])) and row['all_filters_x'] != 1 and row['all_filters_x'] != 2 and row['all_filters_x'] != 3 and (row['ld_x'] or row['ld_x'] is None) else ('ld' if not row['ld_x'] and row['ld_x'] is not None and (row['se_x'] <= .2 or pandas.isnull(row['se_x'])) and (row['se_y'] <= .2 or pandas.isnull(row['se_y'])) and row['all_filters_x'] != 1 and row['all_filters_x'] != 2 and row['all_filters_x'] != 3 else  'FAIL')), row['se_x'] if not pandas.isnull(row['se_x']) else 0, row['se_y'] if not pandas.isnull(row['se_y']) else 0], index = columns), axis = 1)
# else ( else ())
                ndf.columns = ['minicd_p1','minicd_p2', 'minicd_icd1', 'minicd_icd2' , 'minicd_lor1', 'minicd_lor2', 'filter', 'se1', 'se2']
                ndf['ukbb_freq'] = result['ukbb_freq_x']
                ndf['maf'] = result['ukbb_freq_x']
                ndf['HGVSp'] = result['HGVSp_x']
                ndf['HGVSc'] = result['HGVSc_x']
                ndf['rsid'] = result['rsid_x']
                ndf['category'] = result['category_x']
                ndf['variant_identity'] = result['variant_identity']
                ndf['gene_symbol'] = result['gene_symbol_x']
            ndf['shortname1'] = ndf.apply(lambda row: 'NA' if (row['minicd_icd1'] == "NA" or len(str(dfphen[dfphen['title'] == row['minicd_icd1']]['notes'].squeeze()).split(';')) <=  1) else str(dfphen[dfphen['title'] == row['minicd_icd1']]['notes'].squeeze()).split(';')[1].split('=')[1], axis = 1)
            ndf['shortname2'] = ndf.apply(lambda row: 'NA' if (row['minicd_icd2'] == "NA" or len(str(dfphen[dfphen['title'] == row['minicd_icd2']]['notes'].squeeze()).split(';')) <=  1) else str(dfphen[dfphen['title'] == row['minicd_icd2']]['notes'].squeeze()).split(';')[1].split('=')[1], axis = 1)
            ndfmap = ndf.query('filter == "PASS"', inplace = False)
            slope, intercept, r_value, p_value, std_err = stats.linregress(ndfmap['minicd_lor1'],ndfmap['minicd_lor2'])
            ftline = list(slope*ndf['minicd_lor1'] + intercept)
            coverage_stats = ndf.to_dict('records')
            variants_in_gene = coverage_stats
            gbuild = namespace.split('_')[1].lower()
            t = render_template(
                'comparebetas.html',
                variants_in_gene=variants_in_gene,
                coverage_stats=coverage_stats,
                gbuild=gbuild,
                namespace=namespace,
                ftline=ftline,
                p_value=p_value,
                slope=slope,
                cor=r_value)
            return t



@app.route('/<namespace>/gene/interactive/<gene_id>', methods=['GET', 'POST'])
def gene_interactive_page(namespace, gene_id):
    if gene_id in []:
        return open(os.path.join(GENE_CACHE_DIR, '{}.html'.format(gene_id))).read()
    else:
        if request.method == 'POST':
          #  print(dir(request))
          #  print(request.form.keys())
            genes = []
            missense = lof = False
            function = request.form['function']
            if "missense_variant" in function:
                missense = True
            if "lof_variant" in function:
                lof = True
            if request.form['gene_text']:
                text = request.form['gene_text']
                genes = text.split("\n")
        # Remove any empty rows from genes
            for i in range(len(genes)):
                genes[i] = genes[i].rstrip()
                if len(genes[i]) == 0:
                    genes.pop(i)
           # print(genes)
            if 'submit_mrp' in request.form.keys():
                functionfdr = request.form['functionfdr']
                functionphen = request.form.getlist('phenotypes[]')
               # print(functionphen)
                functionfdr = str(functionfdr)
                fdr = functionfdr.strip('[')
                fdr = fdr.strip(']')
                phenidarr = []
                for phenname in functionphen:
                    phenidarr.append(str(phenname))
               # print(phenidarr)
                # Run MRP function
                [key,lbf,clusterval, genedat] = run_mrp(lof=lof, missense=missense, genes=genes, fdr=fdr, phenidarr=phenidarr)
                 # Send results of MRP function to mrp page
                return redirect('/mrp/%s' % key)
            if 'submit_graph' in request.form.keys():
                data = run_mr(lof=lof, missense=missense, genes=genes)
               # print(data['thetainvfindict'])
                key = data['key']
               # print(data.keys())
                graph(**data)
                return redirect('/graph/%s' % key)
            if 'submit_meta' in request.form.keys():
                key = str(random.getrandbits(128))
                annotations = []
                # Add in the selected annotations to the category list                                                                                                                                                                                                                        functionphen = request.form.getlist('phenotypes[]')
                functionphen = request.form['phenotypes']
                phenidarr = []
                phenidarr.append(str(functionphen))
                # From models-example.py
                gene_variant = lookups.get_gene_variant(DB, gene_names=genes, icds=phenidarr)
                keys_M = []
                for gv in gene_variant:
                    keys_M.append('{}-{}-{}-{}'.format(gv['chrom'],
                                                       gv['pos'],
                                                       gv['ref']['val'],
                                                       gv['alt']['val']))
                gene_names_M = gene_variant['gene_name']['val']
                major_consequence_M = gene_variant['consequence']['val']
                HGVSp_M = gene_variant['hgvsp']['val']

                # -- -
                # -- - M x N NumPy Arrays - --
                # -- -
                se_M_N = None
                lor_M_N = None
                for icd in icds:
                    # SciDB lookup
                    # ---
                    variant_icd = lookups.get_variant_icd(DB, gene_names=genes, icds=phenidarr)
                    se = variant_icd['se']['val']
                    if se_M_N is None:
                        se_M_N = se
                    else:
                        se_M_N = numpy.c_[se_M_N, se]
                    lor = variant_icd['lor']['val']
                    if lor_M_N is None:
                        lor_M_N = lor
                    else:
                        lor_M_N = numpy.c_[lor_M_N, lor]
               # print(phenidarr)
                if lof:
                    annotations.append('lof_variant')
                if missense:
                    annotations.append('missense_variant')
                if genes != None:
                    # Find variants with the given annotation
                    b = models.QueryGenome(category=annotations)
                betas, se, pvalues, annotations, protein_annotations, variant_ids, icd, gene_return, rsids, alts, allele_frequencies = b.query_genome(genes,phenidarr)
                betas = [item for sublist in betas for item in sublist]
                se = [item for sublist in se for item in sublist]
                idxdel = []
                for i in range(0,len(se)):
                    if (float(se[i]) == 0) or (float(se[i]) > .5):
                        idxdel.append(i)
                betas = [i for j, i in enumerate(betas) if j not in idxdel]
                protein_annotations = [str(protein_annotations[i].strip('"').strip('[').strip(']').strip("'").split(",")[0].split(":")[1].strip("'")) if len(protein_annotations[i].strip('"').strip('[').strip(']').strip("'").split(",")[0].split(":")) > 1 else 'NA' for i in range(0,len(protein_annotations))]
                variant_ids = [variant_ids[j] + "|" + protein_annotations[j] for j, i in enumerate(variant_ids) if j not in idxdel]
                se = [i for j, i in enumerate(se) if j not in idxdel]
               # meta(key, betas, se, variant_ids,chains = 4, iter = 5000, warmup = 1000, cores = 2)
                return redirect('/meta/%s' % key)
        return get_gene_interactive_page_content(namespace, gene_id)




#@app.route('/<namespace>/geneaggregate/<gene_id>', methods=['GET', 'POST'])
#def geneaggregate_page(namespace, gene_id):
#    return get_geneaggregate_page_content(namespace, gene_id)


@app.route('/<namespace>/coding/geneaggregate/<icd_str>')
def icdaggregate_page(namespace, icd_str):
#    if not check_credentials():
#        return redirect(url_for('login'))
    db = get_db(namespace)
    try:
        # icdlabel = str(icd_str)
        # .strip('ADD').strip('INI').strip('BRMRI')
        # Try several cutoffs to see if there are too many variants to
        # render.  Arbitrary 10k max.
        # passing = False
        cutoff = None
        icd = None

        for p in [.001]:
            icd = lookups.get_icd_variant_by_icd_id_pvalue(db, icd_str, p)
            if len(icd):
                cutoff = p
                # print("CUTOFF",cutoff)
                break
            # print(icd_str,icd)
        print('Rendering ICD10: %s' % icd_str)

        if icd is None or len(icd) == 0:
            icd = [{'Case': 'NA', 'Name': 'NA', 'icd': icd_str}]
        # print(icd_info)

        return render_template(
            'icdgeneaggregate.html',
            icd=icd,
            cutoff=cutoff,
            namespace = namespace
            )
    except Exception as e:
        print('Failed on icd:', icd_str, '; Error=', traceback.format_exc())
        abort(404)



@app.route('/<namespace>/coding/gene-mh/<icd_str>')
def icd_genemh_page(namespace, icd_str):
    db = get_db(namespace)
    try:
        cutoff = None
        icd = None

        for p in [.001]:
            assocset = str(db.list_association_sets()['name'][0])
            df = db.get_phenotype_fields(association_set=assocset, include_pvalue_threshold=True)
            field_identifier = int(df[df['title'] == icd_str]['field_id'])
            pvalthr = float(df[df['title'] == icd_str]['pvalue_threshold'])
            shortname = str(df[df['title'] == icd_str]['notes'].squeeze().split(';')[1].split('=')[1])
            casecnt = str(df[df['title'] == icd_str]['notes'].squeeze().split(';')[0].split('=')[1])
            pvalthr = max(5e-7, pvalthr)
            cuttoff = pvalthr
        icd = [{'Case': casecnt, 'Name': shortname, 'icd': icd_str}]
        # print(icd_info)
        return render_template(
            'gene-mh.html',
            icd=icd,
            cutoff=cutoff,
            namespace=namespace
            )
    except Exception as e:
        print('Failed on icd:', icd_str, '; Error=', traceback.format_exc())
        abort(404)


@app.route('/coding/decomposition-risk/<icd_str>')
def decomposition_risk_page(icd_str):
 #   if not check_credentials():
 #       return redirect(url_for('login'))
 #   elif not check_credentials(permission_group=p_groups.get('rivas-lab')):
 #       abort(404)
    db = get_db()
    try:
        # icdlabel = str(icd_str)
        # .strip('ADD').strip('INI').strip('BRMRI')
        # Try several cutoffs to see if there are too many variants to
        # render.  Arbitrary 10k max.
        # passing = False
        cutoff = None
        icd = None
    except Exception as e:
        print('Failed on decomposition risk page:', icd_str, ';Error=', e)
        abort(404)


@app.route('/<namespace>/transcript/<transcript_id>')
def transcript_page(namespace, transcript_id):
  #  if not check_credentials():
  #      return redirect(url_for('login'))
    db = get_db(namespace)
    try:
        cache_key = 't-transcript-{}'.format(transcript_id)
        t = cache.get(cache_key)
        if t is None:
            transcript = lookups.get_transcript_gene(db, transcript_id)
            transcript_idx = transcript['transcript_idx']
            transcript_id = transcript['transcript_id']
            transcript['exons'] = [
                lookups.cast_pos_info(
                    dict(zip(lookups.EXON_INFO_KEYS, s.split(':'))))
                for s in transcript['exon_info'].split(';')
                if s]

            gene = dict((kg, transcript[kt])
                        for (kt, kg) in (('gene_name', ) * 2,
                                         ('gene_strand', 'strand'),
                                         ('full_gene_name', ) * 2,
                                         ('omim_accession', ) * 2,
                                         ('gene_chrom', 'chrom'),
                                         ('gene_start', 'start'),
                                         ('gene_stop', 'stop'),
                                         ('c_transcript_info', ) * 2,
                                         ('transcript_info', ) * 2,
                                         ('gene_exon_info', 'exon_info'),
                                         ('gene_id', ) * 2))
            gene['transcripts'] = [
                lookups.cast_pos_info(
                    dict(zip(lookups.TRANSCRIPT_INFO_KEYS, s.split(':'))))
                for s in gene['transcript_info'].split(';')
                if s]
            gene['canonical_transcript'] = gene['c_transcript_info'].split(
                ':')[0]

            variants_in_transcript = lookups.get_variants_by_transcript_idx(
                db, transcript_idx, transcript_id)

            coverage_stats = lookups.get_coverage_for_transcript(
                db,
                transcript['xstart'] - EXON_PADDING,
                transcript['xstop'] + EXON_PADDING)

            add_transcript_coordinate_to_variants(
                variants_in_transcript, transcript)

            t = render_template(
                'transcript.html',
                transcript=transcript,
                transcript_json=json.dumps(transcript),
                variants_in_transcript=variants_in_transcript,
                variants_in_transcript_json=json.dumps(variants_in_transcript),
                coverage_stats=coverage_stats,
                coverage_stats_json=json.dumps(coverage_stats),
                gene=gene,
                namespace=namespace,
                gene_json=json.dumps(gene),
            )
            cache.set(cache_key, t, timeout=1000*60)
        print('Rendering transcript: %s' % transcript_id)
        return t
    except Exception as e:
        print('Failed on transcript:', transcript_id, ';Error=', e)
        abort(404)


@app.route('/<namespace>/region/<region_id>')
def region_page(namespace, region_id):
  #  if not check_credentials():
  #      return redirect(url_for('login'))
    db = get_db(namespace)
    try:
        region = region_id.split('-')
        gbuild = namespace.split('_')[1].lower()
        t = None
        if t is None:
            chrom = region[0]
            start = None
            stop = None
            if len(region) == 3:
                chrom, start, stop = region
                start = int(start)
                stop = int(stop)
            if start is None or stop - start > REGION_LIMIT or stop < start:
                return render_template(
                    'region.html',
                    genes_in_region=None,
                    variants_in_region=None,
                    chrom=chrom,
                    start=start,
                    stop=stop,
                    namespace=namespace,
                    gbuild=gbuild
                )
            if start == stop:
                start -= 20
                stop += 20
            genes_in_region = lookups.get_genes_in_region(db, chrom, start, stop).to_dict('records')
            logging.info('genes in region')
            variants_in_region = lookups.get_variants_in_region(db, chrom, start, stop)
            variants_in_region['filter'] = variants_in_region.apply(lambda row: 'PASS' if row['all_filters'] !=  1 and row['all_filters'] != 2 and row['all_filters'] != 3 else "FAIL", axis = 1)
            variants_in_region = variants_in_region.to_dict('records')
            xstart = get_xpos(chrom, start)
            xstop = get_xpos(chrom, stop)
            logging.info('variants in region')
#            coverage_array = lookups.get_coverage_for_bases(db, xstart, xstop)
            t = render_template(
                'region.html',
                genes_in_region=genes_in_region,
                variants_in_region=variants_in_region,
                chrom=chrom,
                start=start,
                stop=stop,
                namespace=namespace,
                gbuild=gbuild
 #               coverage=coverage_array
            )
        logging.info('Rendering region: %s' % region_id)
        return t
    except Exception as e:
        print('Failed on region:', region_id, ';Error=', e)
        abort(404)


@app.route('/<namespace>/dbsnp/<rsid>')
def dbsnp_page(namespace, rsid):
  #  if not check_credentials():
  #      return redirect(url_for('login'))
    db = get_db(namespace)
    try:
        variants = lookups.get_variants_by_rsid(db, rsid)
        chrom = None
        start = None
        stop = None
        print('Rendering rsid: %s' % rsid)
        return render_template(
            'region.html',
            rsid=rsid,
            variants_in_region=variants,
            chrom=chrom,
            start=start,
            stop=stop,
            coverage=None,
            genes_in_region=None,
            namespace=namespace
        )
    except Exception as e:
        print('Failed on rsid:', rsid, ';Error=', e)
        abort(404)



def covid19HGI_read_res_table(static_path='/biobankengine/app/static/', GBE_top=''):
    tbl=os.path.join(static_path, 'covid19HGI/fig/HGIrel5.PRS_PheWAS.1e-5.tsv.gz')
    subset_map={
        'eur_leave_ukbb_23andme': 'EUR (-UKB)',
        'eur_leave_23andme': 'EUR',
        'leave_UKBB_23andme': 'All (-UKB)',
        'leave_23andme': 'All'
    }

    df = pandas.read_csv(tbl, sep='\t')

    df['GBE_short_name'] = ['<a href="{}/RIVAS_HG19/coding/{}">{}</a>'.format(GBE_top, x[0], x[1]) for x in zip(df['GBE_ID'], df['GBE_short_name'])]
    df = df.drop('GBE_ID',axis=1)
    df = df.drop('z_or_t_value',axis=1)
    for col in ['GBE_category']:
        df[col] = df[col].map(lambda x: x.replace('_', ' '))
    for col in ['HGI_suffix']:
        df[col] = df[col].map(lambda x: subset_map[x])
    for col in ['estimate', 'SE']:
        df[col] = df[col].map(lambda x: str(round(x, 3)))
    for col in ['clump_p1']:
        df[col] = df[col].map(lambda x: '{:.0e}'.format(x).replace('-0', '-') )
    for col in ['P']:
        df[col] = df[col].map(lambda x: '{:.2e}'.format(x))

    return(df)

@app.route('/covid19HGI')
def covid19HGI_main():
    df=covid19HGI_read_res_table()

    table_body_str=''.join(['<tr>{}</tr>'.format(
        ''.join(['<td>{}</td>'.format(x) for x in df.iloc[row]])
    ) for row in range(df.shape[0])])
    table_cols=['HGI', 'subset', 'clump p1', 'category', 'trait', 'BETA', 'SE', 'P']
    table_cols_select=['HGI', 'subset', 'clump p1', 'category']

    return render_template(
        'covid19HGI/main.html',
        #=======================#
        _HGI_v = 'HGIrel5',
        _HGI_case_controls=['B2', 'C2', 'B1', 'A2'],
        _HGI_suffices=['eur_leave_ukbb_23andme', 'eur_leave_23andme', 'leave_UKBB_23andme', 'leave_23andme'],
        _clump_p1s=['1e-5', '1e-4', '1e-3'],
        #=======================#
	namespace = 'RIVAS_HG19',
        _HGI_sx = 'eur_leave_ukbb_23andme',
        _HGI_cc = 'B2',
        _clumpp = '1e-5',
        table_cols        = table_cols,
        table_cols_select = table_cols_select,
        table_body_str    = table_body_str
    )

@app.route('/covid19HGI/PRS_PheWAS/<HGI_sx>/')
def covid19HGI_prs_phewas(HGI_sx):
    return render_template(
        'covid19HGI/prs_phewas.html',
        #=======================#
        _HGI_v = 'HGIrel5',
        _HGI_case_controls=['B2', 'C2', 'B1', 'A2'],
        _HGI_suffices=['eur_leave_ukbb_23andme', 'eur_leave_23andme', 'leave_UKBB_23andme', 'leave_23andme'],
        _clump_p1s=['1e-5', '1e-4', '1e-3'],
        #=======================#
	namespace = 'RIVAS_HG19',
        HGI_sx = HGI_sx,
        _HGI_cc = 'B2',
        _clumpp = '1e-5'
    )

@app.route('/covid19HGI/PRS_PheWAS/<HGI_sx>/<HGI_cc>')
def covid19HGI_prs_phewas_more(HGI_sx, HGI_cc):
    return render_template(
        'covid19HGI/prs_phewas_more.html',
        #=======================#
        _HGI_v = 'HGIrel5',
        _HGI_case_controls=['B2', 'C2', 'B1', 'A2'],
        _HGI_suffices=['eur_leave_ukbb_23andme', 'eur_leave_23andme', 'leave_UKBB_23andme', 'leave_23andme'],
        _clump_p1s=['1e-5', '1e-4', '1e-3'],
        _GBE_IDs=['INI30130', 'INI30190', 'INI30610', 'INI10030610'],
        #=======================#
	namespace = 'RIVAS_HG19',
        HGI_sx = HGI_sx,
        HGI_cc = HGI_cc,
        _clumpp = '1e-5'
    )

@app.route('/<namespace>/not_found/<query>')
def not_found_page(namespace, query):
    return render_template(
        'not_found.html',
        query=query,
        namespace=namespace
    )


@app.route('/<namespace>/error/<query>')
@app.errorhandler(404)
def error_page(query):
    namespace='RIVAS_HG19'
    if type(query) == str or type(query) == unicode:
        unsupported = "TTN" if query.upper() in lookups.UNSUPPORTED_QUERIES else None
    else:
        unsupported = None
    return render_template(
        'error.html',
        query=query,
        unsupported=unsupported,
        namespace=namespace
    )


@app.route('/downloads')
def downloads_page():
    namespace='RIVAS_HG19'
    return render_template('downloads.html',namespace=namespace)


@app.route('/about')
def about_page():
    namespace='RIVAS_HG19'
    return render_template('about.html', namespace = namespace)


@app.route('/participants')
def participants_page():
    namespace='RIVAS_HG19'
    return render_template('about.html', namespace=namespace)


@app.route('/terms')
def terms_page():
    namespace = 'RIVAS_HG19'
    return render_template('terms.html', namespace=namespace)


@app.route('/contact')
def contact_page():
    namespace = 'RIVAS_HG19'
    return render_template('contact.html', namespace=namespace)


@app.route('/faq')
def faq_page():
    namespace='RIVAS_HG19'
    return render_template('faq.html', namespace=namespace)


@app.route('/text')
def text_page():
    db = get_db()
    query = request.args.get('text')
    datatype, identifier = lookups.get_awesomebar_result(db, query)
    if datatype in ['gene', 'transcript']:
        gene = lookups.get_gene_by_id(db, identifier)
        link = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr%(chrom)s%%3A%(start)s-%(stop)s" % gene
        output = '''Searched for %s. Found %s.
%s; Canonical: %s.
%s''' % (query, identifier, gene['full_gene_name'], gene['canonical_transcript'], link)
        output += '' if 'omim_accession' not in gene else '''
In OMIM: %(omim_description)s
http://omim.org/entry/%(omim_accession)s''' % gene
        return output
    elif datatype == 'error' or datatype == 'not_found':
        return "Gene/transcript %s not found" % query
    else:
        return "Search types other than gene transcript not yet supported"


@app.route('/gcorr')
def gcorr_page():
#    if not check_credentials():
#        return redirect(url_for('login'))
    namespace='RIVAS_HG19'
    return render_template('gcorr.html',namespace=namespace)

### Dash app ###
# start dash app
dash_app = dash.Dash(__name__, server=app)
dash_app.config.supress_callback_exceptions = True
dash_app.css.append_css({
        "external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"
})

dash_app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    html.Div(id='page-content')
])

# Dash dash_app.callbacks
@dash_app.callback(
    dash.dependencies.Output('pheno-dropdown', 'options'),
    [dash.dependencies.Input('pheno-categories', 'values'),
     dash.dependencies.Input('case-cutoff', 'value'),
    ]
)
def set_possible_phenos(pheno_categories, case_cutoff):
    try:
        case_cutoff = int(case_cutoff)
    except:
        case_cutoff = MIN_CASES
    poss = PHENOS[(PHENOS['GBE_category'].isin(pheno_categories)) &
                  (PHENOS['GBE_N'] >= case_cutoff)]
    return([{'label':x, 'value':x} for x in sorted(poss['phenotype'])])

@dash_app.callback(
    dash.dependencies.Output('gcorr-scatter', 'figure'),
    [dash.dependencies.Input('pheno-dropdown', 'value'),
     dash.dependencies.Input('cluster-method', 'value'),
     dash.dependencies.Input('pheno-categories', 'values'),
     dash.dependencies.Input('z-cutoff', 'value'),
     dash.dependencies.Input('case-cutoff', 'value'),
     dash.dependencies.Input('gcorr-min', 'value'),
     dash.dependencies.Input('gcorr-max', 'value'),
     dash.dependencies.Input('gcorr-radio', 'value'),
     dash.dependencies.Input('pi2-min', 'value'),
     dash.dependencies.Input('pi2-max', 'value'),
     dash.dependencies.Input('pi2-radio', 'value'),
     dash.dependencies.Input('show-zero-estimates', 'values'),
     dash.dependencies.Input('size-var', 'value'),
    ]
)
def update_table(
    selected_phenos,
    cluster_method,
    pheno_categories,
    z_cutoff,
    case_cutoff,
    gcorr_min,
    gcorr_max,
    gcorr_radio,
    pi2_min,
    pi2_max,
    pi2_radio,
    show_zero_estimates,
    size_var,
):
    return(
        gcorr_scatter(
            selected_phenos,
            cluster_method,
            pheno_categories,
            z_cutoff,
            case_cutoff,
            gcorr_min,
            gcorr_max,
            gcorr_radio,
            pi2_min,
            pi2_max,
            pi2_radio,
            show_zero_estimates,
            size_var,
        )
    )

# Dash app callbaack for displaying page. This generally won't be updated.
@dash_app.callback(dash.dependencies.Output('page-content', 'children'),
              [dash.dependencies.Input('url', 'pathname')])
def display_page(pathname):
    if pathname == '/gcorr-app-raw':
        return gcorr_layout


@app.after_request
def apply_caching(response):
    # prevent click-jacking vulnerability identified by BITs
    response.headers["X-Frame-Options"] = "SAMEORIGIN"
    return response


if __name__ == "__main__":
    app.run(host = "0.0.0.0", port = 6000, debug=False, use_debugger=False, use_reloader=True, threaded=True)
