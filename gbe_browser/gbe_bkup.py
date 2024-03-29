from __future__ import division
import sys
import itertools
import collections
import io
import json
import os
import gzip
from parsing import *
import lookups
import random
import models
import numpy
import numpy.matlib
#import rpy2
import pandas
from targeted import mrpmm
from mr import mr
#from graph import graph
from utils import *
from flask import Flask, request, session, g, redirect, url_for, abort, render_template, flash, jsonify
from flask_compress import Compress
from flask_oauth import OAuth
from flask_errormail import mail_on_500
from flask import Response
from collections import defaultdict
from werkzeug.contrib.cache import SimpleCache
from multiprocessing import Process
# STAN workforce
#from optimized import Polygenic
#from optimized import PolygenicCoding
import glob
import traceback
import time
import csv
import logging
from flask_wtf import Form
from wtforms import StringField, BooleanField
from wtforms.validators import DataRequired
import scidbpy
# Plotly dash app packages
from dash_apps import *
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

logging.basicConfig(stream=sys.stderr)
sys.stderr.write("starting gbe.py")
#from pycallgraph import PyCallGraph
#from pycallgraph.output import GraphvizOutput
print(sys.flags.optimize)
def graph(*args):
    pass

def mrp(*args):
    pass

ADMINISTRATORS = (
    'gbe.help@gmail.com',
)

app = Flask(__name__)
app.debug=False

# Get Google API information
# client_secrets.json must be obtained from the Google Developers Console
#with open('client_secrets.json') as secrets:
#    secret_data = json.load(secrets)
#google_client_id = secret_data['web']['client_id']
#google_client_secret = secret_data['web']['client_secret']
#redirect_uri = "/" + secret_data['web']['redirect_uris'][0].split('/')[-1]

#SECRET_KEY = 'development key'
#app.secret_key = SECRET_KEY
#oauth = OAuth()

# Define google oauth protocol
#google = oauth.remote_app('google',
#                          base_url='https://www.google.com/accounts/',
#                          authorize_url='https://accounts.google.com/o/oauth2/auth',
#                          request_token_url=None,
#                          request_token_params={'scope': 'https://www.googleapis.com/auth/userinfo.email',
#                                                'response_type': 'code'},
#                          access_token_url='https://accounts.google.com/o/oauth2/token',
#                          access_token_method='POST',
#                          access_token_params={'grant_type': 'authorization_code'},
#                          consumer_key=google_client_id,
#                          consumer_secret=google_client_secret)


class permission_groups():
    '''
    store dict of set of emails to store white-list of individuals for pages with restricted access.
    To change add/remove/modify the permission groups, please edit permission_groups.tsv
    '''
    def __init__(self, groups=None, input_file=None):
        self._groups = collections.defaultdict(set)
        if groups is not None:
            # one can pass groups as an argument to constructor
            for group, emails in groups:
                self._groups[group].union(set(emails))
        elif input_file is not None:
            # read from file
            with open(input_file) as f:
                # read file and skip empty and comment lines
                input_file_lines = [x for x in f.read().splitlines() if len(x) > 0 and x[0] != '#']
            for line in input_file_lines:
                # the file should be a tsv file. 1st column: group, 2nd column: email
                group, email = line.split('\t')[0], line.split('\t')[1]
                self._groups[group].add(email)
    def get(self, group_name):
        return self._groups[group_name]

p_groups = permission_groups(input_file='permission_groups.tsv')

# Set random seed

mail_on_500(app, ADMINISTRATORS)
Compress(app)
app.config['COMPRESS_DEBUG'] = False
cache = SimpleCache()

GBE_FILES_DIRECTORY = '../gbe_data/'
REGION_LIMIT = 1E5
EXON_PADDING = 50
# Load default config and override config from an environment variable
#test
app.config.update(dict(
    ## Set SCIDB_URL='http://localhost:8080' environment variable
    # SCIDB_URL='http://localhost:8080',
    DEBUG=True,
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

DB = scidbpy.connect(app.config.get('SCIDB_URL', None))
ICD_NAME_MAP = lookups.get_icd_name_map(DB)

# def connect_db():
#     """Connects to the specific database.
#
#     """
#     return scidbpy.connect(app.config.get('SCIDB_URL', None))



def load_icd_info_stats():
    # load ICD10 info
    db = get_db()
    db.icd_info.drop()
    print('Dropped db.icd_info')
    procs = []
    icd_info_file = app.config['ICD_INFO_FILE']
    icd_info = get_icd_info_from_file(icd_info_file)
    db.icd_info.insert(icd_info)
    db.icd_info.ensure_index('icd')
    print('Done Loading ICD INFO!')

def load_gene_models():
    db = get_db()

    db.genes.drop()
    db.transcripts.drop()
    db.exons.drop()
    print('Dropped db.genes, db.transcripts, and db.exons.')

    start_time = time.time()

    canonical_transcripts = {}
    with gzip.open(app.config['CANONICAL_TRANSCRIPT_FILE']) as canonical_transcript_file:
        for gene, transcript in get_canonical_transcripts(canonical_transcript_file):
            canonical_transcripts[gene] = transcript

    omim_annotations = {}
    with gzip.open(app.config['OMIM_FILE']) as omim_file:
        for fields in get_omim_associations(omim_file):
            if fields is None:
                continue
            gene, transcript, accession, description = fields
            omim_annotations[gene] = (accession, description)

    dbnsfp_info = {}
    with gzip.open(app.config['DBNSFP_FILE']) as dbnsfp_file:
        for dbnsfp_gene in get_dbnsfp_info(dbnsfp_file):
            other_names = [other_name.upper() for other_name in dbnsfp_gene['gene_other_names']]
            dbnsfp_info[dbnsfp_gene['ensembl_gene']] = (dbnsfp_gene['gene_full_name'], other_names)

    print('Done loading metadata. Took %s seconds' % int(time.time() - start_time))
    # grab genes from GTF
    start_time = time.time()
    with gzip.open(app.config['GENCODE_GTF']) as gtf_file:
        for gene in get_genes_from_gencode_gtf(gtf_file):
            gene_id = gene['gene_id']
            if gene_id in canonical_transcripts:
                gene['canonical_transcript'] = canonical_transcripts[gene_id]
            if gene_id in omim_annotations:
                gene['omim_accession'] = omim_annotations[gene_id][0]
                gene['omim_description'] = omim_annotations[gene_id][1]
            if gene_id in dbnsfp_info:
                gene['full_gene_name'] = dbnsfp_info[gene_id][0]
                gene['other_names'] = dbnsfp_info[gene_id][1]
                if len(gene['gene_name']) == 0:
                    genesym = ''
                else:
                    genesym = gene['gene_name']
            db.genes.insert(gene, w=0)

    print('Done loading genes. Took %s seconds' % int(time.time() - start_time))

    start_time = time.time()
    db.genes.ensure_index('gene_id')
    db.genes.ensure_index('gene_name_upper')
    db.genes.ensure_index('gene_name')
    db.genes.ensure_index('other_names')
    db.genes.ensure_index('xstart')
    db.genes.ensure_index('xstop')
    print('Done indexing gene table. Took %s seconds' % int(time.time() - start_time))

    # and now transcripts
    start_time = time.time()
    with gzip.open(app.config['GENCODE_GTF']) as gtf_file:
        db.transcripts.insert((transcript for transcript in get_transcripts_from_gencode_gtf(gtf_file)), w=0)
    print('Done loading transcripts. Took %s seconds' % int(time.time() - start_time))

    start_time = time.time()
    db.transcripts.ensure_index('transcript_id')
    db.transcripts.ensure_index('gene_id')
    print('Done indexing transcript table. Took %s seconds' % int(time.time() - start_time))

    # Building up gene definitions
    start_time = time.time()
    with gzip.open(app.config['GENCODE_GTF']) as gtf_file:
        db.exons.insert((exon for exon in get_exons_from_gencode_gtf(gtf_file)), w=0)
    print('Done loading exons. Took %s seconds' % int(time.time() - start_time))

    start_time = time.time()
    db.exons.ensure_index('exon_id')
    db.exons.ensure_index('transcript_id')
    db.exons.ensure_index('gene_id')
    print('Done indexing exon table. Took %s seconds' % int(time.time() - start_time))

    return []



def load_db():
    """
    Load the database
    """
    # Initialize database
    # Don't need to explicitly create tables with mongo, just indices
    confirm = raw_input('This will drop the database and reload. Are you sure you want to continue? [no] ')
    if not confirm.startswith('y'):
        print('Exiting...')
        sys.exit(1)
    all_procs = []
    for load_function in [load_variants_file, load_dbsnp_file, load_base_coverage, load_gene_models, load_icd_stats, load_icd_info_stats]:
        procs = load_function()
        all_procs.extend(procs)
        print("Started %s processes to run %s" % (len(procs), load_function.__name__))
    [p.join() for p in all_procs]
    print('Done! Creating cache...')
    create_cache()
    print('Done!')


def create_cache():
    """
    This is essentially a compile step that generates all cached resources.
    Creates files like autocomplete_entries.txt
    Should be run on every redeploy.
    """
    # create autocomplete_entries.txt
    autocomplete_strings = []
    print >> sys.stderr, "Getting gene names..."
    for gene in get_db().genes.find():
        autocomplete_strings.append(gene['gene_name'])
        if 'other_names' in gene:
            autocomplete_strings.extend(gene['other_names'])
    print >> sys.stderr, "Done! Writing..."
    f = open(os.path.join(os.path.dirname(__file__), 'autocomplete_strings.txt'), 'w')
    for s in sorted(autocomplete_strings):
        f.write(s+'\n')
    f.close()
    print >> sys.stderr, "Done! Getting largest genes..."

    # create static gene pages for genes in
    if not os.path.exists(GENE_CACHE_DIR):
        os.makedirs(GENE_CACHE_DIR)

    # get list of genes ordered by num_variants
    for gene_id in GENES_TO_CACHE:
        try:
            page_content = get_gene_page_content(gene_id)
        except Exception:
            print(Exception)
            continue
        f = open(os.path.join(GENE_CACHE_DIR, '{}.html'.format(gene_id)), 'w')
        f.write(page_content)
        f.close()
    print >> sys.stderr, "Done!"


def precalculate_metrics():
    import numpy
    db = get_db()
    print('Reading %s variants...' % db.variants.count())
    metrics = defaultdict(list)
    binned_metrics = defaultdict(list)
    progress = 0
    start_time = time.time()
    for variant in db.variants.find(fields=['quality_metrics', 'site_quality', 'allele_num', 'allele_count']):
        for metric, value in variant['quality_metrics'].iteritems():
            metrics[metric].append(float(value))
        qual = float(variant['site_quality'])
        metrics['site_quality'].append(qual)
        if variant['allele_num'] == 0: continue
        if variant['allele_count'] == 1:
            binned_metrics['singleton'].append(qual)
        elif variant['allele_count'] == 2:
            binned_metrics['doubleton'].append(qual)
        else:
            for af in AF_BUCKETS:
                if float(variant['allele_count'])/variant['allele_num'] < af:
                    binned_metrics[af].append(qual)
                    break
        progress += 1
        if not progress % 100000:
            print('Read %s variants. Took %s seconds' % (progress, int(time.time() - start_time)))
    print('Done reading variants. Dropping metrics database... ')
    db.metrics.drop()
    print('Dropped metrics database. Calculating metrics...')
    for metric in metrics:
        bin_range = None
        data = map(numpy.log, metrics[metric]) if metric == 'DP' else metrics[metric]
        if metric == 'FS':
            bin_range = (0, 20)
        elif metric == 'VQSLOD':
            bin_range = (-20, 20)
        elif metric == 'InbreedingCoeff':
            bin_range = (0, 1)
        if bin_range is not None:
            data = [x if (x > bin_range[0]) else bin_range[0] for x in data]
            data = [x if (x < bin_range[1]) else bin_range[1] for x in data]
        hist = numpy.histogram(data, bins=40, range=bin_range)
        edges = hist[1]
        # mids = [(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)]
        lefts = [edges[i] for i in range(len(edges)-1)]
        db.metrics.insert({
            'metric': metric,
            'mids': lefts,
            'hist': list(hist[0])
        })
    for metric in binned_metrics:
        hist = numpy.histogram(map(numpy.log, binned_metrics[metric]), bins=40)
        edges = hist[1]
        mids = [(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)]
        db.metrics.insert({
            'metric': 'binned_%s' % metric,
            'mids': mids,
            'hist': list(hist[0])
        })
    db.metrics.ensure_index('metric')
    print('Done pre-calculating metrics!')


def get_db():
    """Opens a new database connection if there is none yet for the
    current application context.

    """
    # if not hasattr(g, 'db_conn'):
    #     g.db_conn = connect_db()
    # return g.db_conn
    return DB

def run_graph():
    key = str(random.getrandbits(128))
    return redirect('/graph/%s' % key)

def run_mrp(lof=True, missense=True, genes=None, fdr=5, phenidarr = ['ICD1462','ICD1463']):
    fdr = int(fdr)/100
    annotations = []

    # Add in the selected annotations to the category list
    if lof:
        annotations.append('lof_variant')
    if missense:
        annotations.append('missense_variant')
    # If the input file has genes
    if genes != None:
        # Find variants with the given annotation
        b = models.QueryGenome(category=annotations)
        key = 't-mrp-' + '_'.join(annotations) + '_' + '_'.join(phenidarr) + '_' + '_'.join(genes)
        # Generate relevant files
# betas, se, pvalues, annotations, protein_annotations, variant_ids, icd, gene_return, rsids, alts, allele_frequencies
        betas, se, pvalues, annotations, protein_annotations, variant_ids, icd, gene_return, rsids, alts, allele_frequencies = b.query_genome(genes,phenidarr)
        t1 = open('test.3.a','w')
        t1.write("BETAS")
        t1.write(betas)
        t1.close()
        numpy.save('s3/' + key + '.betas', betas)
        numpy.save('s3/' + key + '.se', se)
        numpy.save('s3/' + key + '.pvalues', pvalues)
        numpy.save('s3/' + key + '.annotations', annotations)
        numpy.save('s3/' + key + '.protein_annotations', protein_annotations)
        numpy.save('s3/' + key + '.variant_ids', variant_ids)
        numpy.save('s3/' + key + '.icd', icd)
        numpy.save('s3/' + key + '.gene_return', gene_return)
        numpy.save('s3/' + key + '.rsids', rsids)
        numpy.save('s3/' + key + '.alts', alts)
        numpy.save('s3/' + key + '.allele_frequencies', allele_frequencies)
        # Reshape betas from genome query
        C = numpy.matlib.eye(betas.shape[1])
        annotvec = [str(annotations[i].strip('"').strip('[').strip(']').strip("'")) for i in range(0,len(annotations))]

        # Run MRP with 2 clusters, output to the MRP_out subdirectory in the gbe_browser directory
        bicarr = []
        cmax = 4
        fail = 0
        for tmpc in range(1,cmax+1):
            [BIC, AIC, genedat] = mrpmm(betas,se, C, annotvec, gene_return, rsids, variant_ids,tmpc,key,C, numpy.linalg.inv(C),icd, fdr=fdr, niter=151,burn=50,thinning=1,verbose=True, outpath = './MRP_out/')
            if tmpc > 2 and BIC > bicarr[len(bicarr)-1]:
                fail += 1
            if fail >= 2:
                break
            bicarr.append(BIC)

            print(tmpc,BIC,AIC)
#        if tmpc == 5:
#            with PyCallGraph(output=GraphvizOutput()):
#                mrpmm(betas,se, C, annotvec, gene_return, rsids, variant_ids,5,key,C, numpy.linalg.inv(C),icd, fdr=fdr, niter=201,burn=50,thinning=1,verbose=True, outpath = './MRP_out/')
        cminarr = bicarr[1:]
        clustminidx = cminarr.index(min(cminarr))
        clustminval = clustminidx + 2
        clustmaxidx = cminarr.index(max(cminarr))
        clustmaxval = clustmaxidx + 2
        if bicarr[0] > 0 and bicarr[1] > 0:
            lbf2 = 0.5*(bicarr[0] - cminarr[0] )
            lbf = 0.5*( bicarr[0] - cminarr[clustminidx])
        else:
            lbf2 = 0.5*(bicarr[0] - cminarr[0])
            lbf = 0.5*(bicarr[0] - cminarr[clustminidx])
        print("bicstuff",bicarr[0],cminarr[clustminidx])
        lbfout = open("./MRP_out/" + key + ".lbf",'w')
        lbfout.write(str(lbf))
        lbfout.close()
        clustvalue = 1
        if lbf > 1 and (lbf - lbf2) > 1 and lbf2 > 0:
#        if lbf > 1:
            [BIC, AIC, genedat] = mrpmm(betas,se, C, annotvec, gene_return, rsids, variant_ids,clustminval,key,C, numpy.linalg.inv(C), icd, fdr=fdr, niter=1001,burn=500,thinning=1,verbose=True, outpath = './MRP_out/', protectivescan=True)
            print("log bayes factor: ",lbf)
            clustvalue = clustminval
        elif lbf2 > 1:
            clustminval = 2
            [BIC, AIC, genedat] = mrpmm(betas,se, C, annotvec, gene_return, rsids, variant_ids,clustminval,key,C, numpy.linalg.inv(C), icd, fdr=fdr, niter=1001,burn=500,thinning=1,verbose=True, outpath = './MRP_out/', protectivescan=True)
            print("log bayes factor: ",lbf)
            clustvalue = clustminval
            lbf = lbf2
        else:
            [BIC, AIC, genedat] = mrpmm(betas,se, C, annotvec, gene_return, rsids, variant_ids,1,key,C, numpy.linalg.inv(C), icd, fdr=fdr, niter=51,burn=10,thinning=1,verbose=True, outpath = './MRP_out/', protectivescan=True)
        print(bicarr)
        print(icd)
        print("genedat",genedat)
    return([key,lbf,clustvalue, genedat])


def run_mr(lof=True, missense=True, genes=None, phenidarr = ['ICD1462','ICD1463']):
    key = str(random.getrandbits(128))
    annotations = []
    data = {}
    # Add in the selected annotations to the category list
    if lof:
        annotations.append('lof_variant')
    if missense:
        annotations.append('missense_variant')

    # If the input file has genes
    if genes != None:
        # Find variants with the given annotation
        b = models.QueryGenome(category=annotations)

        # Generate relevant files
        betas, se, pvalues, annotations, protein_annotations, variant_ids, icd, gene_return, rsids, alts, allele_frequencies = b.query_genome(genes,phenidarr)

        # Reshape betas from genome query
        C = numpy.matlib.eye(betas.shape[1])
        annotvec = [str(annotations[i].strip('"').strip('[').strip(']').strip("'")) for i in range(0,len(annotations))]
        bicarr = []
        for tmpc in range(1,3):
            returndict = mr(betas,se, C, annotvec, gene_return, rsids, variant_ids,tmpc,key,C, numpy.linalg.inv(C),icd, niter=2001,burn=50,thinning=1,verbose=True, outpath = './MRP_out/')
            BIC = returndict['bic']
            AIC = returndict['aic']
            thetainvdict = returndict['thetainv']
            iter = returndict['iter']
            clustersize = returndict['c']
            bicarr.append(AIC)
        cminarr = bicarr[1:]
        clustminidx = cminarr.index(min(cminarr))
        clustminval = clustminidx + 2
        clustmaxidx = cminarr.index(max(cminarr))
        clustmaxval = clustmaxidx + 2
        if bicarr[0] > 0 and bicarr[1] > 0:
            lbf2 = 0.5*(bicarr[0] - cminarr[0] )
            lbf = 0.5*( bicarr[0] - cminarr[clustminidx])
        else:
            lbf2 = 0.5*(bicarr[0] - cminarr[0])
            lbf = 0.5*(bicarr[0] - cminarr[clustminidx])
       # print("bicstuff",bicarr[0],cminarr[clustminidx])
        lbfout = open("./MRP_out/" + key + ".lbf",'w')
        lbfout.write(str(lbf))
        lbfout.close()
        clustvalue = 1
#        if lbf > 1 and (lbf - lbf2) > 1 and lbf2 > 1:
        if lbf > 1:
            returndict = mr(betas,se, C, annotvec, gene_return, rsids, variant_ids,clustminval,key,C, numpy.linalg.inv(C), icd, niter=2001,burn=1000,thinning=1,verbose=True, outpath = './MRP_out/')
            BIC = returndict['bic']
            AIC = returndict['aic']
            thetainvdict = returndict['thetainv']
            iter = returndict['iter']
            clustersize = returndict['c']
          #  print("log bayes factor: ",lbf)
            clustvalue = clustersize
            thetainvfindict = thetainvdict
           # print(returndict['thetainv'], "HERE")
        elif lbf2 > 1:
            clustminval = 2
            returndict = mr(betas,se, C, annotvec, gene_return, rsids, variant_ids,clustminval,key,C, numpy.linalg.inv(C), icd, niter=2001,burn=1000,thinning=1,verbose=True, outpath = './MRP_out/')
            BIC = returndict['bic']
            AIC = returndict['aic']
            thetainvdict = returndict['thetainv']
            iter = returndict['iter']
            clustersize = returndict['c']
           # print("log bayes factor: ",lbf)
            clustvalue = clustersize
            thetainvfindict = thetainvdict
           # print(returndict['thetainv'], "HERE2")
        else:
            returndict = mr(betas,se, C, annotvec, gene_return, rsids, variant_ids,2,key,C, numpy.linalg.inv(C), icd, niter=2001,burn=1000,thinning=1,verbose=True, outpath = './MRP_out/')
            BIC = returndict['bic']
            AIC = returndict['aic']
            thetainvdict = returndict['thetainv']
            iter = returndict['iter']
            clustersize = returndict['c']
            thetainvfindict = thetainvdict
            clustvalue = clustersize
            thetainvfindict = thetainvdict
           # print(returndict['thetainv'], "HERE3")
        data = {}
        data['key'] = key
        data['lbf'] = lbf
        data['thetainvfindict'] = numpy.array(thetainvfindict).tolist()
        data['clustvalue'] = clustvalue
        data['iter'] = iter
      #  print(bicarr)
      #  print(icd)
      #  print(data['thetainvfindict'])
    return data


# @app.teardown_appcontext
# def close_db(error):
#     """Closes the database again at the end of the request."""
#     if hasattr(g, 'db_conn'):
#         g.db_conn.close()

@app.route('/')
def homepage():
#    if check_credentials():
#        return redirect(url_for('search_page'))
#    return render_template('homepage.html')
    return render_template('search_page.html')

@app.route('/search')
def search_page():
#    if not check_credentials():
#        return redirect(url_for('login'))
    return render_template('search_page.html')


@app.route('/autocomplete/<query>')
def awesome_autocomplete(query):
    if not hasattr(g, 'autocomplete_strings'):
        g.autocomplete_strings = [s.strip() for s in open(os.path.join(os.path.dirname(__file__), 'autocomplete_strings.txt'))]
    suggestions = lookups.get_awesomebar_suggestions(g, query)
    return Response(json.dumps([{'value': s} for s in suggestions]),  mimetype='application/json')


@app.route('/awesome')
def awesome():
    db = get_db()
    query = request.args.get('query')
    datatype, identifier = lookups.get_awesomebar_result(db, query)

    print("Searched for %s: %s" % (datatype, identifier))
    if datatype == 'gene':
        return redirect('/gene/{}'.format(identifier))
    elif datatype == 'transcript':
        return redirect('/transcript/{}'.format(identifier))
    elif datatype == 'dbsnp_variant_set':
        return redirect('/dbsnp/{}'.format(identifier))
    elif datatype == 'variant':
        return redirect('/variant/{}'.format(identifier))
    elif datatype == 'icd10':
        return redirect('/coding/{}'.format(identifier))
    elif datatype == 'region':
        return redirect('/region/{}'.format(identifier))
    elif datatype == 'error':
        return redirect('/error/{}'.format(identifier))
    elif datatype == 'not_found':
        return redirect('/not_found/{}'.format(identifier))
    else:
        raise Exception


@app.route('/variant/<variant_str>')
def variant_icd_page(variant_str):
  #  if not check_credentials():
  #      return redirect(url_for('login'))
    db = get_db()
    try:
        chrom, pos = variant_str.split('-')
        # pos, ref, alt = get_minimal_representation(pos, ref, alt)
        variant = lookups.get_variant_ann_by_chrom_pos(db, chrom, pos)
        if variant is None:
            variant = {
                'chrom': chrom,
                'pos': pos,
                'xpos': get_xpos(chrom, int(pos))
            }
        print(variant.keys())
        consequences = None
        ordered_csqs = None
        if 'vep_annotations' in variant:
            variant['vep_annotations'] = order_vep_by_csq(variant['vep_annotations'])  # Adds major_consequence
            ordered_csqs = [x['major_consequence'] for x in variant['vep_annotations']]
            ordered_csqs = reduce(lambda x, y: ','.join([x, y]) if y not in x else x, ordered_csqs, '').split(',') # Close but not quite there
            consequences = defaultdict(lambda: defaultdict(list))
            for annotation in variant['vep_annotations']:
                annotation['HGVS'] = get_proper_hgvs(annotation)
                consequences[annotation['major_consequence']][annotation['Gene']].append(annotation)
       # print(xpos)
        icdstats = lookups.get_icd_by_chrom_pos(db, chrom, pos)
        indexes = []
        seend = {}
        for idx in range(0,len(icdstats)):
            # ICD10=T81/Name=Complications of procedures, not elsewhere classified/Chapter=T/L95OR=0.97/U95OR=2.04/OR=1.40/pvalue=0.0756463/l10pval=1.12/Case=1229
            item = icdstats[idx]
            icd10 = item['icd']
            item['Code'] = icd10
            # icd10info = lookups.get_icd_info(db, icd10)
            if 'Name' not in item:
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
                # item['Name'] = icd10info[0]['Name']
                item['Group'] = icd10[0] # default value
                groups = ['RH', 'FH', 'HC', 'cancer', 'ADD', 'INI', 'MED', 'BIN', 'BRMRI', 'BROADBIN', 'BROADQT', 'INI_FC', 'BIN_FC']
                for group in groups:
                    if icd10.startswith(group):
                        item['Group'] = group
                        break
                item['OR'] = format(float(item['or_val']), '.4g')
                item['LOR'] = format(float(item['lor']), '.4g')
                item['L95OR'] = format(float(item['l95or']), '.4g')
                item['U95OR'] = format(float(item['u95or']), '.4g')
                item['pvalue'] = format(float(item['pvalue']), '.4g')
                item['l10pval'] = format(float(item['log10pvalue']), '.4g')
                item['SE'] = format(float(item['se']), '.4g')
                if float(item['pvalue']) == 0:
                    item['pvalue'] = numpy.finfo(float).eps
                    item['pvalue'] = format(float(item['pvalue']),'.4g')
                    item['l10pval'] = 250
                # item['Case'] = icd10info[0]['Case']
                se =  format(float(item['se']), '.4g')
                if float(item['l10pval']) < 1 or float(se) >= .5 or (float(se) >= .08 and item['OR'] == item['LOR']) or int(item['Case']) <= 100  or item['Code'] == "HC67" or icd10 in seend:
                    indexes.append(idx)
                seend[icd10] = icd10
        for index in sorted(indexes, reverse=True):
            del icdstats[index]
        print('Rendering variant: %s' % variant_str)
        return render_template(
            'variant.html',
            variant=variant,
            icdstats=icdstats,
            consequences=consequences,
            ordered_csqs=ordered_csqs,
            debug_message='',
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
        if group in set(['INI', 'INI_FC', 'BROADQT']):
            beta_or_lor = 'BETA'
            beta_or_lor_val = plot_d_dict[group]['LOR']
            beta_or_lor_l95 = np.array(plot_d_dict[group]['LOR']) - 1.96 * np.array(plot_d_dict[group]['SE'])
            beta_or_lor_u95 = np.array(plot_d_dict[group]['LOR']) + 1.96 * np.array(plot_d_dict[group]['SE'])

        else:
            beta_or_lor = 'OR'
            beta_or_lor_val = np.exp(np.array(plot_d_dict[group]['LOR']))
            beta_or_lor_l95 = np.exp(np.array(plot_d_dict[group]['LOR']) - 1.96 * np.array(plot_d_dict[group]['SE']))
            beta_or_lor_u95 = np.exp(np.array(plot_d_dict[group]['LOR']) + 1.96 * np.array(plot_d_dict[group]['SE']))

        group_len = len(plot_d_dict[group]['icd'])
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
        'x':    plot_d_dict[group]['icd'],
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
        'x':    plot_d_dict[group]['icd'],
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


@app.route('/coding/<icd_str>')
def icd_page(icd_str):
#    if not check_credentials():
#        return redirect(url_for('login'))
    db = get_db()
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
            'icd.html',
            icd=icd,
            cutoff=cutoff
            )
    except Exception as e:
        print('Failed on icd:', icd_str, '; Error=', traceback.format_exc())
        abort(404)


@app.route('/codinguncertain/<icd_str>')
def icduncertain_page(icd_str):
  #  if not check_credentials():
  #      return redirect(url_for('login'))
    db = get_db()
    try:
        # icdlabel = str(icd_str)
        # .strip('ADD').strip('INI').strip('BRMRI')
        # Try several cutoffs to see if there are too many variants to
        # render.  Arbitrary 10k max.
        # passing = False
        cutoff = None
        icd = None

        for p in [.001]:
            icd = lookups.get_icd_variant_by_icd_id_pvalue_uncertain(db, icd_str, p)
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
            'icduncertain.html',
            icd=icd,
            cutoff=cutoff
            )
    except Exception as e:
        print('Failed on icd:', icd_str, '; Error=', traceback.format_exc())
        abort(404)



@app.route('/coding/gene/<icd_str>')
def codinggene_page(icd_str):
  #  if not check_credentials():
  #      return redirect(url_for('login'))
    db = get_db()
    try:
        # icdlabel = str(icd_str)
        # .strip('ADD').strip('INI').strip('BRMRI')
        # Try several cutoffs to see if there are too many variants to
        # render.  Arbitrary 10k max.
        # passing = False
        cutoff = None
        icd = None
        for p in [.00001]:
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
            'icdgene.html',
            icd=icd,
            cutoff=cutoff
            )
    except Exception as e:
        print('Failed on icd:', icd_str, '; Error=', traceback.format_exc())
        abort(404)



@app.route('/coding/phenotype/<icd_str>')
def codingphenotype_page(icd_str):
  #  if not check_credentials():
  #      return redirect(url_for('login'))
    db = get_db()
    try:
        # icdlabel = str(icd_str)
        # .strip('ADD').strip('INI').strip('BRMRI')
        # Try several cutoffs to see if there are too many variants to
        # render.  Arbitrary 10k max.
        # passing = False
        cutoff = None
        icd = None
        for p in [.00001]:
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
            'phenotype.html',
            icd=icd,
            cutoff=cutoff
            )
    except Exception as e:
        print('Failed on icd:', icd_str, '; Error=', traceback.format_exc())
        abort(404)


@app.route('/target/1')
def target_page():
  #  if not check_credentials():
  #      return redirect(url_for('login'))
    db = get_db()
    try:
        passing = False
        cutoff = None

        for p in [.0000001]:
            icd = lookups.get_icd_variant_by_pvalue(db, p)
            print("icd here",icd)
            cutoff = p
#            if len(icd):
#                break

        #     vars = lookups.get_significant_prot(db, p, 0)
        #     print(vars)
        #     if len(vars) < 10000000000:
        #         passing = True
        #     if passing:
        #         cutoff = p
        #         break
        # variants = get_prot_variant_table(0,cutoff)

        return render_template(
            'target.html',
            icd = icd,
            cutoff=cutoff
            )
    except Exception as e:
        print('Failed on protective scan Error=', traceback.format_exc())
        abort(404)


@app.route('/power')
def power_page():
  #  if not check_credentials():
  #      return redirect(url_for('login'))
    try:
        return render_template(
            'power.html'
            )
    except Exception as e:
        print('Failed on protective scan Error=', traceback.format_exc())
        abort(404)


#@app.route('/decomposition-dev')
@app.route('/degas')
def decomposition_dev_page():    
  #  if not check_credentials():
  #      return redirect(url_for('login'))
  #  elif not check_credentials(permission_group=p_groups.get('rivas-lab')):
  #      abort(404)

    db = get_db()
    
    try:
        # this list file will be automatically updated by cron job that calls the following script:
        # https://github.com/rivas-lab/decomposition/blob/master/src/decomposition_dataset_list.sh
        with open('./static/decomposition/decomposition_datasets.lst') as f:
            dataset_list = f.read().splitlines()
                
        return render_template(
            'decomposition-dev.html',    
            dataset_list = dataset_list,
        )
    
    except Exception as e:
        print('Unknown Error=', traceback.format_exc())
        abort(404)


@app.route('/decomposition-internal/<dataset>')
def decomposition_internal_page(dataset):    
  #  if not check_credentials():
  #      return redirect(url_for('login'))
  #  elif not check_credentials(permission_group=p_groups.get('rivas-lab')):
  #      abort(404)
    db = get_db()
    
    init_idx_pc  = 0
    init_idx_phe = 0
    init_idx_var = 0
    debug_str = 'debug'
    
    try:
        return render_template(
            'decomposition-internal.html',    
            init_idx_pc  = init_idx_pc,
            init_idx_phe = init_idx_phe,
            init_idx_var = init_idx_var,
            dataset = dataset,
            debug_str = debug_str
        )            
    
    except Exception as e:
        print('Unknown Error=', traceback.format_exc())
        abort(404)


@app.route('/decomposition-app')
def decomposition_app_page():    
    abort(404)

def decomposition_app_page_null():
 #   if not check_credentials():
 #       return redirect(url_for('login'))
    db = get_db()
    
    try:
        with open('./static/decomposition/decomposition-app.lst') as f:
            dataset_list = f.read().splitlines()
                
        return render_template(
            'decomposition-app.html',    
            dataset_list = dataset_list,
        )
    
    except Exception as e:
        print('Unknown Error=', traceback.format_exc())
        abort(404)
        
        

@app.route('/decomposition/<dataset>')
def decomposition_page(dataset):    
 #   if not check_credentials():
 #       return redirect(url_for('login'))
 #   elif not check_credentials(permission_group=p_groups.get('rivas-lab')):
 #       abort(404)
    db = get_db()
    
    
    if(dataset == 'PTVs'):
        init_idx_pc  = 6
        init_idx_phe = 429
        init_idx_var = 8503    
    else:
        init_idx_pc  = 0
        init_idx_phe = 0
        init_idx_var = 0
        
    debug_str = 'debug'
    
    try:
        if(dataset == "20170930_EMBL-Stanford_coding-nonMHC_z"):
            return render_template(
                'decomposition-20170930.html',    
                init_idx_pc  = init_idx_pc,
                init_idx_phe = init_idx_phe,
                init_idx_var = init_idx_var,
                dataset = dataset,
                debug_str = debug_str
            )
        elif(dataset[:8] == "20171011"):
            return render_template(
                'decomposition-20171011.html',    
                init_idx_pc  = init_idx_pc,
                init_idx_phe = init_idx_phe,
                init_idx_var = init_idx_var,
                dataset = dataset,
                debug_str = debug_str
            )
        else:
            return render_template(
                'decomposition.html',    
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
 #   if not check_credentials():
 #       return redirect(url_for('login'))
    db = get_db()

    try:
        return render_template('hla-assoc.html')

    except Exception as e:
        print('Unknown Error=', traceback.format_exc())
        abort(404)


def get_icd_variant_table(icd, p):
    db = get_db()
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


def get_prot_variant_table(lor, p):
    db = get_db()
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



@app.route('/intensity/<affy_str>')
def intensity_page(affy_str):
 #   if not check_credentials():
 #       return redirect(url_for('login'))
    db = get_db()
    try:
        n_UKBL = 11
        n_UKBB = 95
        
        
        print('Rendering Intensity page: %s' % affy_str)
        variant = lookups.get_icd_affyid(db, affy_str)
        # print(variant)
        return render_template(
            'intensity.html',
            affy=affy_str,
            variant=variant,
            UKBL_idx = [1 + x for x in range(n_UKBL)],
            UKBB_idx = [1 + x for x in range(n_UKBB)]
            )
    except Exception as e:
        print('Failed on affy id:', affy_str, '; Error=', traceback.format_exc())
        abort(404)

@app.route('/gene/<gene_id>', methods=['GET', 'POST'])
def gene_page(gene_id):
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
        return get_gene_page_content(gene_id)


@app.route('/gene/interactive/<gene_id>', methods=['GET', 'POST'])
def gene_interactive_page(gene_id):
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
        return get_gene_interactive_page_content(gene_id)




@app.route('/geneaggregate/<gene_id>', methods=['GET', 'POST'])
def geneaggregate_page(gene_id):
    return get_geneaggregate_page_content(gene_id)


@app.route('/coding/geneaggregate/<icd_str>')
def icdaggregate_page(icd_str):
#    if not check_credentials():
#        return redirect(url_for('login'))
    db = get_db()
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
            cutoff=cutoff
            )
    except Exception as e:
        print('Failed on icd:', icd_str, '; Error=', traceback.format_exc())
        abort(404)



@app.route('/coding/gene-mh/<icd_str>')
def icd_genemh_page(icd_str):
 #   if not check_credentials():
 #       return redirect(url_for('login'))
    db = get_db()
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
            'gene-mh.html',
            icd=icd,
            cutoff=cutoff
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
            'decomposition-risk.html',
            icd=icd,
            cutoff=cutoff
            )
    except Exception as e:
        print('Failed on icd:', icd_str, '; Error=', traceback.format_exc())
        abort(404)







@app.route('/runPolyCoding', methods = ['GET', 'POST'])
def runPolyCoding_page():
    if request.method == 'POST':
        db = get_db()
        genes = []
        missense = lof = False
        function = request.form['function']
        if "missense_variant" in function:
            missense = True
        if "lof_variant" in function:
            lof = True
        annotations = []
        # Add in the selected annotations to the category list
        if lof:
            annotations.append('lof_variant')
        if missense:
            annotations.append('missense_variant')
        print(annotations)
        functionphen = request.form.getlist('phenotypes[]')
        print(functionphen)
        phenidarr = []
        for phenname in functionphen:
            phenidarr.append(str(phenname))
        # Run Poly function
        if 'submit_polygenic' in request.form.keys():
            functionphen = request.form.getlist('phenotypes[]')
            phenidarr = []
            for phenname in functionphen:
                        phenidarr.append(str(phenname))
            phenidarr.sort()
            keycheck = '_'.join(phenidarr)
            key = keycheck
            if os.path.exists("static/images/PolygenicCoding/PolygenicCoding_" + str(keycheck) + ".svg"):
                return redirect("/polygeniccoding/%s" % key)
            else:
                b = models.QueryGenome(category=annotations)
                # Generate relevant files
                betas, se, pvalues, annotations, protein_annotations, variant_ids, icd, gene_return, rsids, alts, allele_frequencies = b.query_genome(None, phenidarr)
                labels = phenidarr
                print(betas)
                criteria = ~(se==0).any(1)
                print(criteria)
                betas = betas[criteria]
                se = se[criteria]
                print(labels)
                print(key)
                print(betas)
                print(len(betas))
                print(betas.shape)
                print(se)
                print(len(se))
                f5test = open('test5.tmp','w')
                f5test.write("BETAS")
                f5test.write(str(len(betas)))
                f5test.write(key)
                f5test.close()
               # PolygenicCoding(key, betas, se, labels, chains = 8, iter = 200, warmup = 100, cores = 8)
        return redirect('/polygeniccoding/%s' % key)
    form = LoginForm()
    if form.validate_on_submit():
        flash('Login requested for OpenID="%s", remember_me=%s' %
              (form.openid.data, str(form.remember_me.data)))
        return redirect('/runPolyCoding.html')
    return render_template('runPolyCoding.html')

@app.route('/runMRP', methods=['GET', 'POST'])
def runMRP_page():

    if request.method == 'POST':

        genes = []

        missense = lof = False
        function = request.form['function']
        functionfdr = request.form['functionfdr']
        functionfdr = str(functionfdr)
        functionphen = request.form.getlist('phenotypes[]')
        phenidarr = []
        for phenname in functionphen:
            phenidarr.append(str(phenname))
        fdr = functionfdr.strip('[')
        fdr = fdr.strip(']')
        if "missense_variant" in function:
            missense = True
        if "lof_variant" in function:
            lof = True

        if request.form['gene_text']:
            text = request.form['gene_text']
            genes = text.split("\n")


        '''
        # Files being a huge pain
        if 'file' not in request.files:
            flash('No file part')
            return str(request.files)
            return redirect(request.url)
        file = request.files['file']
        if file.filename == '':
            return "file blank"
            flash('No selected file')
            return redirect(request.url)
        if file and allowed_file(file.filename):
            return "Found file"
            filename = secure_filename(file.filename)
            new_file = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            file.save(new_file)
            stuff = []
            with open(new_file) as f:
                for line in f:
                    stuff.append(f.rstrip())
            return(",".join(stuff))

            #genes, lof=True
            #results = run_mrp(genes, lof, etc)
            return redirect('/mrp.html')

            #cs = ComputeScore(twenty_three_file=new_file, score_file=scoring_matrix, score_file_index=scoring_index)
            #r = cs.twenty_three_and_me_analysis()
            #return results(r)
        '''

        # Remove any empty rows from genes
        for i in range(len(genes)):
            genes[i] = genes[i].rstrip()
            if len(genes[i]) == 0:
                genes.pop(i)
        # Caching the page
        annotations = []
        # Add in the selected annotations to the category list
        if lof:
            annotations.append('lof_variant')
        if missense:
            annotations.append('missense_variant')
        key = '_'.join(annotations) + '_' + '_'.join(phenidarr) + '_' + '_'.join(genes)
        cache_key = 't-mrp-{}'.format(key)
        t = cache.get(cache_key)
        # if not cached
        if t is None:
            if os.path.exists(os.path.join('MRP_cache','{}.html'.format(cache_key))):
                return  open(os.path.join('MRP_cache','{}.html'.format(cache_key))).read()
            elif os.path.exists(os.path.join("MRP_out/",key,".mcmc.posteriors")):
                pass
            else:
                [key,lbf,clusterval, genedat] = run_mrp(lof=lof, missense=missense, genes=genes, fdr=fdr, phenidarr=phenidarr)
        print("REndering mrp page: %s" % cache_key)        # Send results of MRP function to mrp page
        return redirect('/mrp/%s' % cache_key)


    form = LoginForm()
    if form.validate_on_submit():
        flash('Login requested for OpenID="%s", remember_me=%s' %
              (form.openid.data, str(form.remember_me.data)))
        return redirect('/runMRP.html')
    return render_template('runMRP.html')


def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

class LoginForm(Form):
    openid = StringField('openid', validators=[DataRequired()])
    remember_me = BooleanField('remember_me', default=False)
    age = StringField('age', validators=[DataRequired()])

def get_gene_page_content(gene_id):
 #   if not check_credentials():
 #       return redirect(url_for('login'))
    db = get_db()
    try:
        cache_key = 't-gene-{}'.format(gene_id)
        t = cache.get(cache_key)
        if t is None:
            gene = lookups.get_gene_by_id(db, gene_id)
            if gene is None:
                abort(404)
            #print(gene_id)
            gene_idx = gene['gene_idx']

            variants_in_gene = lookups.get_variants_by_gene_idx(db, gene_idx, gene_id)

            transcripts_in_gene = [
                lookups.cast_pos_info(
                    dict(zip(lookups.TRANSCRIPT_INFO_KEYS, s.split(':'))))
                for s in gene['transcript_info'].split(';')
                if s]

            #print(variants_in_gene)
            # Get some canonical transcript and corresponding info

            transcript_id = gene['canonical_transcript']
            transcript_idx = gene['transcript_idx']
            transcript = lookups.add_xpos(
                lookups.cast_pos_info(
                    dict(zip(lookups.TRANSCRIPT_INFO_KEYS,
                             gene['c_transcript_info'].split(':')))))
            transcript['transcript_id'] = transcript_id
            transcript['exons'] = [
                lookups.cast_pos_info(
                    dict(zip(lookups.EXON_INFO_KEYS, s.split(':'))))
                for s in gene['exon_info'].split(';')
                if s]

            variants_in_transcript = lookups.get_variants_by_transcript_idx(
                db, transcript_idx, transcript_id)

            coverage_stats = lookups.get_coverage_for_transcript(
                db,
                transcript['xstart'] - EXON_PADDING,
                transcript['xstop'] + EXON_PADDING)

            add_transcript_coordinate_to_variants(
                variants_in_transcript, transcript)

            #print("\n\n\n\nVariants in gene")
            #print(variants_in_gene[0])

            # Add minicd info
            for variant in variants_in_gene:
                variant["minicd_info"] = ICD_NAME_MAP.get(variant["minicd"],
                                                          "info_not_found")

            #print ("Transcripts in gene")
            #print(transcripts_in_gene)

            #print ("Variants in transcript")
            #print (variants_in_transcript)

            t = render_template(
                'gene.html',
                gene=gene,
                transcript=transcript,
                variants_in_gene=variants_in_gene,
                variants_in_transcript=variants_in_transcript,
                transcripts_in_gene=transcripts_in_gene,
                coverage_stats=coverage_stats
            )
            cache.set(cache_key, t, timeout=1000*60)
        print('Rendering gene: %s' % gene_id)
        return t
    except Exception as e:
        print('Failed on gene:', gene_id, ';Error=', e)
        abort(404)


def get_gene_interactive_page_content(gene_id):
 #   if not check_credentials():
 #       return redirect(url_for('login'))
    db = get_db()
    try:
        cache_key = 't-gene-{}'.format(gene_id)
        t = None
        if t is None:
            gene = lookups.get_gene_by_id(db, gene_id)
            if gene is None:
                abort(404)
            #print(gene_id)
            gene_idx = gene['gene_idx']

            variants_in_gene = lookups.get_variants_by_gene_idx(db, gene_idx, gene_id)

            transcripts_in_gene = [
                lookups.cast_pos_info(
                    dict(zip(lookups.TRANSCRIPT_INFO_KEYS, s.split(':'))))
                for s in gene['transcript_info'].split(';')
                if s]

            #print(variants_in_gene)
            # Get some canonical transcript and corresponding info

            transcript_id = gene['canonical_transcript']
            transcript_idx = gene['transcript_idx']
            transcript = lookups.add_xpos(
                lookups.cast_pos_info(
                    dict(zip(lookups.TRANSCRIPT_INFO_KEYS,
                             gene['c_transcript_info'].split(':')))))
            transcript['transcript_id'] = transcript_id
            transcript['exons'] = [
                lookups.cast_pos_info(
                    dict(zip(lookups.EXON_INFO_KEYS, s.split(':'))))
                for s in gene['exon_info'].split(';')
                if s]

            variants_in_transcript = lookups.get_variants_by_transcript_idx(
                db, transcript_idx, transcript_id)

            coverage_stats = lookups.get_coverage_for_transcript(
                db,
                transcript['xstart'] - EXON_PADDING,
                transcript['xstop'] + EXON_PADDING)

            add_transcript_coordinate_to_variants(
                variants_in_transcript, transcript)

            #print("\n\n\n\nVariants in gene")
            #print(variants_in_gene[0])

            # Add minicd info
            for variant in variants_in_gene:
                variant["minicd_info"] = ICD_NAME_MAP.get(variant["minicd"],
                                                          "info_not_found")

            #print ("Transcripts in gene")
            #print(transcripts_in_gene)

            #print ("Variants in transcript")
            #print (variants_in_transcript)

            t = render_template(
                'geneinteractive.html',
                gene=gene,
                transcript=transcript,
                variants_in_gene=variants_in_gene,
                variants_in_transcript=variants_in_transcript,
                transcripts_in_gene=transcripts_in_gene,
                coverage_stats=coverage_stats
            )
            cache.set(cache_key, t, timeout=1000*60)
        print('Rendering gene: %s' % gene_id)
        return t
    except Exception as e:
        print('Failed on gene:', gene_id, ';Error=', e)
        abort(404)




def get_geneaggregate_page_content(gene_id):
  #  if not check_credentials():
  #      return redirect(url_for('login'))
    db = get_db()
    try:
        t = None
        if t is None:
            gene = lookups.get_gene_by_id(db, gene_id)
            if gene is None:
                abort(404)
            #print(gene_id)
            gene_idx = gene['gene_idx']

            variants_in_gene = lookups.get_variants_by_gene_idx(db, gene_idx, gene_id)

            transcripts_in_gene = [
                lookups.cast_pos_info(
                    dict(zip(lookups.TRANSCRIPT_INFO_KEYS, s.split(':'))))
                for s in gene['transcript_info'].split(';')
                if s]

            #print(variants_in_gene)
            # Get some canonical transcript and corresponding info

            transcript_id = gene['canonical_transcript']
            transcript_idx = gene['transcript_idx']
            transcript = lookups.add_xpos(
                lookups.cast_pos_info(
                    dict(zip(lookups.TRANSCRIPT_INFO_KEYS,
                             gene['c_transcript_info'].split(':')))))
            transcript['transcript_id'] = transcript_id
            transcript['exons'] = [
                lookups.cast_pos_info(
                    dict(zip(lookups.EXON_INFO_KEYS, s.split(':'))))
                for s in gene['exon_info'].split(';')
                if s]

            variants_in_transcript = lookups.get_variants_by_transcript_idx(
                db, transcript_idx, transcript_id)

            coverage_stats = lookups.get_coverage_for_transcript(
                db,
                transcript['xstart'] - EXON_PADDING,
                transcript['xstop'] + EXON_PADDING)

            add_transcript_coordinate_to_variants(
                variants_in_transcript, transcript)

            #print("\n\n\n\nVariants in gene")
            #print(variants_in_gene[0])

            # Add minicd info
            for variant in variants_in_gene:
                variant["minicd_info"] = ICD_NAME_MAP.get(variant["minicd"],
                                                          "info_not_found")

            #print ("Transcripts in gene")
            #print(transcripts_in_gene)

            #print ("Variants in transcript")
            #print (variants_in_transcript)

            t = render_template(
                'geneaggregate.html',
                gene=gene,
                transcript=transcript,
                variants_in_gene=variants_in_gene,
                variants_in_transcript=variants_in_transcript,
                transcripts_in_gene=transcripts_in_gene,
                coverage_stats=coverage_stats
            )
        print('Rendering gene: %s' % gene_id)
        return t
    except Exception as e:
        print('Failed on gene:', gene_id, ';Error=', e)
        abort(404)





@app.route('/transcript/<transcript_id>')
def transcript_page(transcript_id):
  #  if not check_credentials():
  #      return redirect(url_for('login'))
    db = get_db()
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
                gene_json=json.dumps(gene),
            )
            cache.set(cache_key, t, timeout=1000*60)
        print('Rendering transcript: %s' % transcript_id)
        return t
    except Exception as e:
        print('Failed on transcript:', transcript_id, ';Error=', e)
        abort(404)


@app.route('/region/<region_id>')
def region_page(region_id):
  #  if not check_credentials():
  #      return redirect(url_for('login'))
    db = get_db()
    try:
        region = region_id.split('-')
        cache_key = 't-region-{}'.format(region_id)
        t = cache.get(cache_key)
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
                    coverage=None
                )
            if start == stop:
                start -= 20
                stop += 20
            genes_in_region = lookups.get_genes_in_region(db, chrom, start, stop)
            variants_in_region = lookups.get_variants_in_region(db, chrom, start, stop)
            xstart = get_xpos(chrom, start)
            xstop = get_xpos(chrom, stop)
            coverage_array = lookups.get_coverage_for_bases(db, xstart, xstop)
            t = render_template(
                'region.html',
                genes_in_region=genes_in_region,
                variants_in_region=variants_in_region,
                chrom=chrom,
                start=start,
                stop=stop,
                coverage=coverage_array
            )
        print('Rendering region: %s' % region_id)
        return t
    except Exception as e:
        print('Failed on region:', region_id, ';Error=', e)
        abort(404)


@app.route('/dbsnp/<rsid>')
def dbsnp_page(rsid):
  #  if not check_credentials():
  #      return redirect(url_for('login'))
    db = get_db()
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
            genes_in_region=None
        )
    except Exception as e:
        print('Failed on rsid:', rsid, ';Error=', e)
        abort(404)


@app.route('/not_found/<query>')
def not_found_page(query):
    return render_template(
        'not_found.html',
        query=query
    )


@app.route('/error/<query>')
@app.errorhandler(404)
def error_page(query):
    if type(query) == str or type(query) == unicode:
        unsupported = "TTN" if query.upper() in lookups.UNSUPPORTED_QUERIES else None
    else:
        unsupported = None
    return render_template(
        'error.html',
        query=query,
        unsupported=unsupported
    )


@app.route('/downloads')
def downloads_page():
    return render_template('downloads.html')


@app.route('/about')
def about_page():
    return render_template('about.html')


@app.route('/participants')
def participants_page():
    return render_template('about.html')


@app.route('/terms')
def terms_page():
    return render_template('terms.html')


@app.route('/contact')
def contact_page():
    return render_template('contact.html')


@app.route('/faq')
def faq_page():
    return render_template('faq.html')


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

@app.route('/graph/<key>')
def graphpage(key):
    if not check_credentials():
        return redirect(url_for('login'))
    db = get_db()
    try:
        return render_template(
            'graph.html',
            graph_key=key
            )
    except Exception as e:
        print('Failed: %s' % e)
        abort(404)

@app.route('/polygenic/<key>')
def polygenicpage(key):
    if not check_credentials():
        return redirect(url_for('login'))
    db = get_db()
    try:
        return render_template(
            'polygenic.html',
            polygenic_key=key
            )
    except Exception as e:
        print('Failed: %s' % e)
        abort(404)

@app.route('/polygeniccoding/<key>')
def polygeniccodingpage(key):
    if not check_credentials():
        return redirect(url_for('login'))
    db = get_db()
    try:
        return render_template(
            'polygeniccoding.html',
            polygenic_key=key
            )
    except Exception as e:
        print('Failed: %s' % e)
        abort(404)

@app.route('/meta/<key>')
def metapage(key):
    if not check_credentials():
        return redirect(url_for('login'))
    db = get_db()
    try:
        return render_template(
            'meta.html',
            meta_key=key
            )
    except Exception as e:
        print('Failed: %s' % e)
        abort(404)

@app.route('/mrp/<key>')
def mrp(key):
    if not check_credentials():
        return redirect(url_for('login'))
    t = cache.get(key)
    if t is not None:
        return t
    db = get_db()
    with open("./MRP_out/"+key+".mcmc.posteriors", "r") as inFile:
        probabilities = []
        admixture_data = []
        variants = []
        nullmembership = []
        for line in inFile.readlines():
            line = line.strip()
            var_info = line.split("\t")
            varid = var_info[0]
            var = var_info[4]
            variant = {}
            varidpage = varid.split('-')[0] + '-' + varid.split('-')[1]
            variant["variant"] = varidpage
            variant["varid"] = var
            for j in range(5, len(var_info)):
                variant["%s" % (j-4)] = var_info[j]
            if float(var_info[5]) >= .99:
                        continue
            nullmembership.append(float(var_info[5]))
            variants.append(variant)
        idxnewarr = [b[0] for b in sorted(enumerate(nullmembership),key=lambda i:i[1], reverse = True)]
        variants = [variants[idxnewarr[i]] for i in range(0,len(idxnewarr))]
        admixture_data.append(variants)
    with open("./MRP_out/"+key+".mcmc.gene.posteriors", "r") as inFile:
        admixture_datagene = []
        genes = []
        for line in inFile.readlines():
            line = line.strip()
            gene_info = line.split("\t")
            geneid = gene_info[0]
            gene = {}
            gene["gene"] = geneid
            for j in range(1, int((len(gene_info)-1)/3)+1):
                gene["%s" % (j-1)] = gene_info[j]
            genes.append(gene)
        admixture_datagene.append(genes)
    with open("./MRP_out/" + key + ".lbf","r") as inFile:
        for line in inFile.readlines():
            line = line.strip()
            line = line.split()
            l10bf = "{0:.3g}".format(float(line[0])/numpy.log(10))
    with open("./MRP_out/"+key+".mcmc.bc", "r") as inFile:
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
            line = line.strip()
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
                    tmp.append(float(clust_info[j]))
                if j % 3 == 2:
                    tmp.append(float(clust_info[j]))
                if j % 3 == 0:
                    tmp.append(float(clust_info[j]))
                    cluster.append(tmp)
            cluster_data.append([cluster_num, cluster])
    fdrf = open("./MRP_out/"+key+".fdr", "r").readlines()
    fdrnum = fdrf[0].rstrip().split()[0]
    fdr_data = []
    for line in fdrf[1:]:
        line = line.strip()
        line = line.split()
        fdr_data.append(line[0])
    try:
        t = render_template(
            'mrp.html',
            plot_data=admixture_data,
            gene_data=admixture_datagene,
            num_figs=cluster_data,
            log10_bf=l10bf,
            cluster_num=str(int(cluster_num)+1),
            fdr=fdrnum,
            fdr_data=fdr_data
   #         phenid_arr=headerarr
            )
        cache.set(key, t, timeout = 0)
        f = io.open(os.path.join('MRP_cache','{}.html'.format(str(key))), 'w', encoding="utf-8")
        f.write(t)
        f.close()
        return t
    except Exception as e:
        print('Failed: %s' % e)
        abort(404)

@app.route('/gcorr')
def gcorr_page():
#    if not check_credentials():
#        return redirect(url_for('login'))
    return render_template('gcorr.html')

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
    poss = PHENOS[(PHENOS['category'].isin(pheno_categories)) & 
                  (PHENOS['numcases'] >= case_cutoff)]
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

@app.route('/login')
def login():
    callback = url_for('authorized', _external=True)
    return google.authorize(callback=callback)

@app.route('/login_page')
def login_page():
    return render_template('login_page.html')


#@app.route(redirect_uri)
#@google.authorized_handler
def authorized(resp):
    access_token = resp['access_token']
    session['access_token'] = access_token, ''
    return redirect(url_for('search_page'))

#@google.tokengetter
def get_acess_token():
    return session.get('access_token')

# Function to check Google Login Credentials for current session
def check_credentials(permission_group=None):
    # Debugging.  Sick of logging in and getting redirected
    if app.debug:
       return True
    access_token = session.get('access_token')
    if access_token is None:
        return False

    access_token = access_token[0]
    from urllib2 import Request, urlopen, URLError

    headers = {'Authorization': 'OAuth ' + access_token}
    req = Request('https://www.googleapis.com/oauth2/v1/userinfo',
                  None, headers)
    try:
        res = urlopen(req)
        data = json.loads(res.read())
        # todo log email address data['email']
        if ((permission_group is None) or (data['email'] in permission_group)):
            return True
        else:
            return False
        #return True
    except URLError as e:
        if e.code == 401:
            # Unauthorized - bad token
            session.pop('access_token', None)
            return redirect(url_for('login_page'))
        return False


@app.after_request
def apply_caching(response):
    # prevent click-jacking vulnerability identified by BITs
    response.headers["X-Frame-Options"] = "SAMEORIGIN"
    return response

if __name__ == "__main__":
    app.run(host = "0.0.0.0", port = 5000, debug=False, use_debugger=False, use_reloader=False)
