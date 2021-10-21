from __future__ import division
import sys
import itertools
import json
import os
import pymongo
import pysam
import gzip
from parsing import *
import lookups_mongodb
import random
#import models_mongodb
from utils_mongodb import *
from flask import Flask, request, session, g, redirect, url_for, abort, render_template, flash, jsonify
from flask_compress import Compress
from flask import Flask, redirect, url_for, session
#from flask_oauth import OAuth
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
import pandas
import numpy
import numpy.matlib
import numpy as np
from flask_wtf import Form
from wtforms import StringField, BooleanField
from wtforms.validators import DataRequired
logging.basicConfig(stream=sys.stderr)
sys.stderr.write("starting gbe.py")
#from pycallgraph import PyCallGraph
#from pycallgraph.output import GraphvizOutput
print(sys.flags.optimize)
def graph(*args):
    pass

def mrp(*args):
    pass

import collections
from collections import defaultdict
#from dash_apps import *
#import dash
#import dash_core_components as dcc
#import dash_html_components as html
#import plotly.graph_objs as go

ADMINISTRATORS = (
    'ukbb.browsers.errors@gmail.com',
)

app = Flask(__name__)
app.debug=True

SECRET_KEY = 'development key'
app.secret_key = SECRET_KEY



# Set random seed


mail_on_500(app, ADMINISTRATORS)
Compress(app)
app.config['COMPRESS_DEBUG'] = True
cache = SimpleCache()

GBE_FILES_DIRECTORY = '../gbe_data/'
REGION_LIMIT = 1E5
EXON_PADDING = 50
# Load default config and override config from an environment variable
#test
app.config.update(dict(
    DB_HOST='localhost',
    DB_PORT=27017, 
    DB_NAME='gbe', 
    DEBUG=True,
    SECRET_KEY='development key',
    LOAD_DB_PARALLEL_PROCESSES = 8,  # contigs assigned to threads, so good to make this a factor of 24 (eg. 2,3,4,6,8)
    SITES_VCFS=glob.glob(os.path.join(os.path.dirname(__file__), GBE_FILES_DIRECTORY, 'icd10ukbb.ukbiobank.merge.sort.vcf.gz')),
    GENCODE_GTF=os.path.join(os.path.dirname(__file__), GBE_FILES_DIRECTORY, 'gencode.gtf.gz'),
    CANONICAL_TRANSCRIPT_FILE=os.path.join(os.path.dirname(__file__), GBE_FILES_DIRECTORY, 'canonical_transcripts.txt.gz'),
    OMIM_FILE=os.path.join(os.path.dirname(__file__), GBE_FILES_DIRECTORY, 'omim_info.txt.gz'),
    BASE_COVERAGE_FILES=glob.glob(os.path.join(os.path.dirname(__file__), GBE_FILES_DIRECTORY, 'coverage', 'Panel2016.*.coverage.txt.gz')),
    BASE_ICDSTATS_FILES=glob.glob(os.path.join(os.path.dirname(__file__), GBE_FILES_DIRECTORY, 'icdstats', 'Panel*.icdstats.txt.gz')),
    ICD_INFO_FILE=os.path.join(os.path.dirname(__file__), GBE_FILES_DIRECTORY, 'icdstats', 'icdinfo.txt'),
    ICD_STATS_FILES=glob.glob(os.path.join(os.path.dirname(__file__), GBE_FILES_DIRECTORY, 'icdassoc','hybrid','*c*.hybrid.rewritewna.gz')),
    QT_STATS_FILES=glob.glob(os.path.join(os.path.dirname(__file__), GBE_FILES_DIRECTORY, 'icdassoc','hybrid','*.linear.rewrite.gz')),
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

def connect_db():
    """
    Connects to the specific database.
    """
    print(app.config['DB_HOST'])
    print(app.config['DB_PORT'])
    client = pymongo.MongoClient(host=app.config['DB_HOST'], port=app.config['DB_PORT'])
    print(client)
    return client[app.config['DB_NAME']]


def parse_tabix_file_subset(tabix_filenames, subset_i, subset_n, record_parser):
    """
    Returns a generator of parsed record objects (as returned by record_parser) for the i'th out n subset of records
    across all the given tabix_file(s). The records are split by files and contigs within files, with 1/n of all contigs
    from all files being assigned to this the i'th subset.

    Args:
        tabix_filenames: a list of one or more tabix-indexed files. These will be opened using pysam.Tabixfile
        subset_i: zero-based number
        subset_n: total number of subsets
        record_parser: a function that takes a file-like object and returns a generator of parsed records
    """
    start_time = time.time()
   # open_tabix_files = []
    open_tabix_files = [pysam.Tabixfile(tabix_filename) for tabix_filename in tabix_filenames]
    tabix_file_contig_pairs = [(tabix_file, contig) for tabix_file in open_tabix_files for contig in tabix_file.contigs]
    tabix_file_contig_subset = tabix_file_contig_pairs[subset_i : : subset_n]  # get every n'th tabix_file/contig pair
    short_filenames = ", ".join(map(os.path.basename, tabix_filenames))
    num_file_contig_pairs = len(tabix_file_contig_subset)
    print(("Loading subset %(subset_i)s of %(subset_n)s total: %(num_file_contig_pairs)s contigs from "
           "%(short_filenames)s") % locals())
    counter = 0
    for tabix_file, contig in tabix_file_contig_subset:
        header_iterator = tabix_file.header
        records_iterator = tabix_file.fetch(contig, 0, 10**9)
        for parsed_record in record_parser(itertools.chain(header_iterator, records_iterator)):
            counter += 1
            yield parsed_record

            if counter % 100000 == 0:
                seconds_elapsed = int(time.time()-start_time)
                print(("Loaded %(counter)s records from subset %(subset_i)s of %(subset_n)s from %(short_filenames)s "
                       "(%(seconds_elapsed)s seconds)") % locals())

    print("Finished loading subset %(subset_i)s from  %(short_filenames)s (%(counter)s records)" % locals())


def parse_tabix_file_subset2(tabix_filenames, subset_i, subset_n, record_parser):
    """
    Returns a generator of parsed record objects (as returned by record_parser) for the i'th out n subset of records
    across all the given tabix_file(s). The records are split by files and contigs within files, with 1/n of all contigs
    from all files being assigned to this the i'th subset.

    Args:
        tabix_filenames: a list of one or more tabix-indexed files. These will be opened using pysam.Tabixfile
        subset_i: zero-based number
        subset_n: total number of subsets
        record_parser: a function that takes a file-like object and returns a generator of parsed records
    """
    start_time = time.time()
#    open_tabix_files = []
    open_tabix_files = [pysam.Tabixfile(tabix_filename) for tabix_filename in tabix_filenames]
    tabix_file_contig_pairs = [(tabix_file, contig) for tabix_file in open_tabix_files for contig in tabix_file.contigs]
    tabix_file_contig_subset = tabix_file_contig_pairs[subset_i : : subset_n]  # get every n'th tabix_file/contig pair
    short_filenames = ", ".join(map(os.path.basename, tabix_filenames))
    num_file_contig_pairs = len(tabix_file_contig_subset)
    print(("Loading subset %(subset_i)s of %(subset_n)s total: %(num_file_contig_pairs)s contigs from "
           "%(short_filenames)s") % locals())
    counter = 0
    for tabix_file, contig in tabix_file_contig_subset:
        header_iterator = tabix_file.header
        records_iterator = tabix_file.fetch(contig, 0, 10**9)
        for parsed_record in record_parser(itertools.chain(header_iterator, records_iterator), tabix_filenames):
            counter += 1
            yield parsed_record

            if counter % 100000 == 0:
                seconds_elapsed = int(time.time()-start_time)
                print(("Loaded %(counter)s records from subset %(subset_i)s of %(subset_n)s from %(short_filenames)s "
                       "(%(seconds_elapsed)s seconds)") % locals())
        tabix_file.close()
    print("Finished loading subset %(subset_i)s from  %(short_filenames)s (%(counter)s records)" % locals())


def load_base_coverage():
    def load_coverage(coverage_files, i, n, db):
        coverage_generator = parse_tabix_file_subset(coverage_files, i, n, get_base_coverage_from_file)
        try:
            db.base_coverage.insert(coverage_generator, w=0)
        except pymongo.errors.InvalidOperation:
            pass  # handle error when coverage_generator is empty

    db = get_db()
    db.base_coverage.drop()
    print("Dropped db.base_coverage")
    # load coverage first; variant info will depend on coverage
    db.base_coverage.ensure_index('xpos')

    procs = []
    coverage_files = app.config['BASE_COVERAGE_FILES']
    num_procs = app.config['LOAD_DB_PARALLEL_PROCESSES']
    random.shuffle(app.config['BASE_COVERAGE_FILES'])
    for i in range(num_procs):
        p = Process(target=load_coverage, args=(coverage_files, i, num_procs, db))
        p.start()
        procs.append(p)
    return procs

def load_base_icdstats():
    def load_icdstats(icdstats_files, i, n, db):
        # done
        icdstats_generator = parse_tabix_file_subset(icdstats_files, i, n, get_base_icdstats_from_file)
        try:
            db.base_icdstats.insert(icdstats_generator, w=0)
        except pymongo.errors.InvalidOperation:
            pass  # handle error when icdstats_generator is empty

    db = get_db()
    # update this 
    db.base_icdstats.drop()
    print("Dropped db.base_icdstats")
    # load coverage first; variant info will depend on coverage
    db.base_icdstats.ensure_index('xpos')

    procs = []
    icdstats_files = app.config['BASE_ICDSTATS_FILES']
    num_procs = app.config['LOAD_DB_PARALLEL_PROCESSES']
    random.shuffle(app.config['BASE_ICDSTATS_FILES'])
    for i in range(num_procs):
        p = Process(target=load_icdstats, args=(icdstats_files, i, num_procs, db))
        p.start()
        procs.append(p)
    return procs

def load_variants_file():
    def load_variants(sites_file, i, n, db):
        variants_generator = parse_tabix_file_subset([sites_file], i, n, get_variants_from_sites_vcf)
        try:
            db.variants.insert(variants_generator, w=0)
        except pymongo.errors.InvalidOperation:
            pass  # handle error when variant_generator is empty

    db = get_db()
    db.variants.drop()
    print("Dropped db.variants")

    # grab variants from sites VCF
    db.variants.ensure_index('xpos')
    db.variants.ensure_index('xstart')
    db.variants.ensure_index('xstop')
    db.variants.ensure_index('rsid')
    db.variants.ensure_index('genes')
    db.variants.ensure_index('transcripts')

    sites_vcfs = app.config['SITES_VCFS']
    if len(sites_vcfs) > 1:
        raise Exception("More than one sites vcf file found: %s" % sites_vcfs)

    procs = []
    num_procs = app.config['LOAD_DB_PARALLEL_PROCESSES']
    for i in range(num_procs):
        p = Process(target=load_variants, args=(sites_vcfs[0], i, num_procs, db))
        p.start()
        procs.append(p)
    return procs

    #print 'Done loading variants. Took %s seconds' % int(time.time() - start_time)

def load_icd_stats():
    def load_icd(icd_files, i, n, db):
        icd_generator = parse_tabix_file_subset2(icd_files, i, n, get_icd_from_file)
        # add get_icd_from_file
        try:
            db.icd.insert(icd_generator, w=0)
        except pymongo.errors.InvalidOperation:
            pass
    def load_qt(qt_files, i, n, db):
        qt_generator = parse_tabix_file_subset2(qt_files, i, n, get_qt_from_file)
        # add get_qt_from_file
        try:
            db.icd.insert(qt_generator, w=0)
        except pymongo.errors.InvalidOperation:
            pass
    # load ICD10 code
    # allows us to search for ICDK50 , ICDK51 in the search bar -> then open icd/K50 presents a manhattan plot - allow points to link to variant page
    db = get_db()
    db.icd.drop()
    print('Dropped db.icd')
    db.icd.ensure_index('icd')
    db.icd.ensure_index('xpos')
    db.icd.ensure_index('icdind')
    db.icd.ensure_index('affyid')
    procs = []
    icd_files = app.config['ICD_STATS_FILES']
    qt_files = app.config['QT_STATS_FILES']
    num_procs = app.config['LOAD_DB_PARALLEL_PROCESSES']
    print(icd_files)
    keyarr = []
    qtkeyarr = []
    keylist = {}
    qtkeylist = {}
    for filecd in icd_files:
        keyarr.append(os.path.basename(filecd).split('.')[1])
        icdk = os.path.basename(filecd).split('.')[1]
        if icdk in keylist:
            keylist[icdk].append(filecd)
        else:
            keylist[icdk] = [filecd]
    for filecd in qt_files:
        qtkeyarr.append(os.path.basename(filecd).split('.')[1])
        icdk = os.path.basename(filecd).split('.')[1]
        if icdk in qtkeylist:
            qtkeylist[icdk].append(filecd)
        else:
            qtkeylist[icdk] = [filecd]
    keyarr = list(set(keyarr))
    qtkeyarr = list(set(qtkeyarr))
    #random.shuffle(app.config['ICD_STATS_FILES'])
    
    for icdkey in keyarr:
        for i in range(num_procs):
            p = Process(target=load_icd, args=(keylist[icdkey], i, num_procs, db))
            p.start()
            procs.append(p)
        [p.join() for p in procs]
    procs = []
    for icdkey in qtkeyarr:
        for i in range(num_procs):
            p = Process(target=load_qt, args=(qtkeylist[icdkey], i, num_procs, db))
            p.start()
            procs.append(p)
        [p.join() for p in procs]
    return []


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

    gwasdict = {}
    lofdict = {}
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
    print(transcript for transcript in get_transcripts_from_gencode_gtf(gtf_file))
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


def load_dbsnp_file():
    db = get_db()

    def load_dbsnp(dbsnp_file, i, n, db):
        if os.path.isfile(dbsnp_file + ".tbi"):
            dbsnp_record_generator = parse_tabix_file_subset([dbsnp_file], i, n, get_snp_from_dbsnp_file)
            try:
                db.dbsnp.insert(dbsnp_record_generator, w=0)
            except pymongo.errors.InvalidOperation:
                pass  # handle error when coverage_generator is empty

        else:
            with gzip.open(dbsnp_file) as f:
                db.dbsnp.insert((snp for snp in get_snp_from_dbsnp_file(f)), w=0)

    db.dbsnp.drop()
    db.dbsnp.ensure_index('rsid')
    db.dbsnp.ensure_index('xpos')
    start_time = time.time()
    dbsnp_file = app.config['DBSNP_FILE']

    print("Loading dbsnp from %s" % dbsnp_file)
    if os.path.isfile(dbsnp_file + ".tbi"):
        num_procs = app.config['LOAD_DB_PARALLEL_PROCESSES']
    else:
        # see if non-tabixed .gz version exists
        if os.path.isfile(dbsnp_file):
            print(("WARNING: %(dbsnp_file)s.tbi index file not found. Will use single thread to load dbsnp."
                "To create a tabix-indexed dbsnp file based on UCSC dbsnp, do: \n"
                "   wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp141.txt.gz \n"
                "   gzcat snp141.txt.gz | cut -f 1-5 | bgzip -c > snp141.txt.bgz \n"
                "   tabix -0 -s 2 -b 3 -e 4 snp141.txt.bgz") % locals())
            num_procs = 1
        else:
            raise Exception("dbsnp file %s(dbsnp_file)s not found." % locals())

    procs = []
    for i in range(num_procs):
        p = Process(target=load_dbsnp, args=(dbsnp_file, i, num_procs, db))
        p.start()
        procs.append(p)

    return procs
    #print 'Done loading dbSNP. Took %s seconds' % int(time.time() - start_time)

    #start_time = time.time()
    #db.dbsnp.ensure_index('rsid')
    #print 'Done indexing dbSNP table. Took %s seconds' % int(time.time() - start_time)


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
    """
    Opens a new database connection if there is none yet for the
    current application context.
v    """
    if not hasattr(g, 'db_conn'):
        g.db_conn = connect_db()
    return g.db_conn

def run_graph():
    key = str(random.getrandbits(128))
    return redirect('/graph/%s' % key)


# @app.teardown_appcontext
# def close_db(error):
#     """Closes the database again at the end of the request."""
#     if hasattr(g, 'db_conn'):
#         g.db_conn.close()

@app.route('/')
def homepage():
    return render_template('search_page.html')

@app.route('/search')
def search_page():
    return render_template('search_page.html')


@app.route('/autocomplete/<query>')
def awesome_autocomplete(query):
    if not hasattr(g, 'autocomplete_strings'):
        g.autocomplete_strings = [s.strip() for s in open(os.path.join(os.path.dirname(__file__), 'autocomplete_strings.txt'))]
    suggestions = lookups_mongodb.get_awesomebar_suggestions(g, query)
    return Response(json.dumps([{'value': s} for s in suggestions]),  mimetype='application/json')


@app.route('/awesome')
def awesome():
    db = get_db()
    query = request.args.get('query')
    datatype, identifier = lookups_mongodb.get_awesomebar_result(db, query)

    print("Searched for %s: %s" % (datatype, identifier))
    if datatype == 'gene':
        return redirect('/gene/{}'.format(identifier))
    elif datatype == 'transcript':
        return redirect('/transcript/{}'.format(identifier))
    elif datatype == 'variant':
        return redirect('/variant/{}'.format(identifier))
    elif datatype == 'icd10':
        return redirect('/coding/{}'.format(identifier))
    elif datatype == 'region':
        return redirect('/region/{}'.format(identifier))
    elif datatype == 'dbsnp_variant_set':
        return redirect('/dbsnp/{}'.format(identifier))
    elif datatype == 'error':
        return redirect('/error/{}'.format(identifier))
    elif datatype == 'not_found':
        return redirect('/not_found/{}'.format(identifier))
    else:
        raise Exception


@app.route('/variant/<variant_str>')
def variant_icd_page(variant_str):
    db = get_db()
    try:
        chrom, pos, ref, alt = variant_str.split('-')
        pos = int(pos)
        # pos, ref, alt = get_minimal_representation(pos, ref, alt)
        xpos = get_xpos(chrom, pos)
        variant = lookups_mongodb.get_variant(db, xpos)
        if variant is None:
            variant = {
                'chrom': chrom,
                'pos': pos,
                'xpos': xpos
            }
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
        icdstats = lookups_mongodb.get_variant_icd(db, xpos)
        indexes = []
        for idx in range(0,len(icdstats)): 
            # ICD10=T81/Name=Complications of procedures, not elsewhere classified/Chapter=T/L95OR=0.97/U95OR=2.04/OR=1.40/pvalue=0.0756463/l10pval=1.12/Case=1229
            item = icdstats[idx]
            icd10 = item['icd']
            item['Code'] = icd10
            icd10info = lookups_mongodb.get_icd_info(db, icd10)
            print(icd10info)
            print("ICD10 INFO!!")
            if len(icd10info) == 0:
                item['Name'] = 'NA'
                item['Group'] = 'NA'
                item['OR'] = 1
                item['LOR'] = 1
                item['L95OR'] = 1
                item['SE'] = 1
                item['U95OR'] = 1
                item['pvalue'] = 1
                item['l10pval'] = 0
                item['log10pvalue'] = 0
                item['Case'] = 'NA'
                indexes.append(idx)
            else:
                item['Name'] = icd10info[0]['Name']
                if icd10[0:2] == "RH":
                    item['Group'] = "RH"
                elif icd10[0:2] == "FH":
                    item['Group'] = "FH"
                elif icd10[0:2] == "HC":
                    item['Group'] = "HC"
                elif icd10[0:6] == "cancer":
                    item['Group'] = "cancer"
                elif icd10[0:3] == "ADD":
                    item['Group'] = "ADD"
                elif icd10[0:3] == "INI":
                    item['Group'] = "INI"
                elif icd10[0:5] == "BRMRI":
                    item['Group'] = "BRMRI"
                else:
                    item['Group'] = icd10[0]
                item['OR'] = format(float(item['stats'][0]['or']), '.4g')
                item['LOR'] = format(float(item['stats'][0]['lor']), '.4g')
                item['L95OR'] = format(float(item['stats'][0]['l95or']), '.4g')
                item['U95OR'] = format(float(item['stats'][0]['u95or']), '.4g')
                item['SE'] = format(float(item['stats'][0]['se']), '.4g')
                item['pvalue'] = format(float(item['stats'][0]['pvalue']), '.4g')
                item['l10pval'] = format(float(item['stats'][0]['log10pvalue']), '.4g')
                item['Case'] = icd10info[0]['Case']
                item['log10pvalue'] = item['l10pval']
                se =  format(float(item['stats'][0]['se']), '.4g') 
                if float(item['l10pval']) <= 1 or float(se) >= .5:
                    indexes.append(idx)
        #for index in sorted(indexes, reverse=True):
        #    del icdstats[index]
        print('MAF')
        variant['maf'] = variant['allele_freq']
        print('Rendering variant: %s' % variant_str)
        print("VARIANT STRING")
        print(icdstats)
        return render_template(
            'variant_mongodb.html',
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
        arrtmp = np.array([x[key] for x in icdstats])
	if len(arrtmp.shape) == 2 and arrtmp.shape[1] == 1:
            print(arrtmp)
            arrtmp = np.squeeze(arrtmp)
        plot_d_raw[key] = arrtmp
    print(pandas.DataFrame(plot_d_raw))
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
        if group not in set(['Disease_outcome', 'HC','RH', 'FH', 'cancer', 'Mental_health', 'Family_history','Others', 'Health_and_medical_history', 'Psychosocial_factors', 'Digestive_health']):
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

@app.route('/coding/<icd_str>')
def icd_page(icd_str):
    db = get_db()
    try:
        icdlabel = icd_str
        #.strip('ADD').strip('INI').strip('BRMRI')
        # Try several cutoffs to see if there are too many variants to render.  Arbitrary 10k max.
        passing = False
        cutoff = None
        icd = None

        for p in [.001, .0001, .00001]:
            icd = lookups_mongodb.get_icd_significant(db, str(icd_str), p)
            if len(icd) < 100000:
                passing = True
            if passing:
                cutoff = p
               # print("CUTOFF",cutoff)
                break
        icd_info = lookups_mongodb.get_icd_info(db, str(icdlabel))
        variants = get_icd_variant_table(icd_str, cutoff)
        print('Rendering ICD10: %s' % icdlabel)
        if len(icd_info) == 0:
            icd_info.append({'Case': 'NA', 'Name': 'NA', 'icd': icdlabel})

        return render_template(
            'icd_mongodb.html',
            icd=icd,
            variants_in_gene = variants,
            icd_info=icd_info,
            cutoff=cutoff
            )
    except Exception as e:
        print('Failed on icd:', icd_str, '; Error=', traceback.format_exc())
        abort(404)


def get_icd_variant_table(icd, p):
    db = get_db()
    significant_variant_ids = {}
    significant_variants = lookups_mongodb.get_icd_significant(db, icd, p)
    for v in significant_variants:
        significant_variant_ids[v['xpos']] = {}
        significant_variant_ids[v['xpos']]['pvalue'] = v['stats'][0]['pvalue']
        significant_variant_ids[v['xpos']]['or'] = v['stats'][0]['or']
        significant_variant_ids[v['xpos']]['log10pvalue'] = v['stats'][0]['log10pvalue']
        significant_variant_ids[v['xpos']]['lor_val'] = v['stats'][0]['lor']
    variants = lookups_mongodb.get_variants_by_id(db, significant_variant_ids.keys())
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
        v['log10pvalue'] = significant_variant_ids[v['xpos']]['log10pvalue']
        v['lor_val'] = significant_variant_ids[v['xpos']]['lor_val']
    return variants


@app.route('/intensity/<affy_str>')
def intensity_page(affy_str):
    db = get_db()
    try:
        print('Rendering Intensity page: %s' % affy_str)
        variant = lookups_mongodb.get_variant_affy(db, affy_str)
       # print(variant)
        return render_template(
            'intensity.html',
            affy=affy_str, 
            variant=variant
            )
    except Exception as e:
        print('Failed on affy id:', affy_str, '; Error=', traceback.format_exc())
        abort(404)

@app.route('/gene/<gene_id>', methods=['GET', 'POST'])
def gene_page(gene_id):
    if gene_id in GENES_TO_CACHE:
        print(gene_id)
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
        return get_gene_page_content(gene_id)





def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

class LoginForm(Form):
    openid = StringField('openid', validators=[DataRequired()])
    remember_me = BooleanField('remember_me', default=False)
    age = StringField('age', validators=[DataRequired()])

def get_gene_page_content(gene_id):
    db = get_db()
    try:
        gene = lookups_mongodb.get_gene_by_name(db, gene_id)
        print(gene)
        if gene is None:
            abort(404)
        cache_key = 't-gene-{}'.format(gene_id)
        t = cache.get(cache_key)
        if t is None:
           # print(gene_id)
            variants_in_gene = lookups_mongodb.get_variants_in_region(db, gene['chrom'], gene['start'] - EXON_PADDING, gene['stop'] + EXON_PADDING)
            transcripts_in_gene = lookups_mongodb.get_transcripts_in_gene(db, gene_id)
            # Get some canonical transcript and corresponding info
            print(variants_in_gene)
            transcript_id = gene['canonical_transcript']
            print(transcript_id)
            transcript = lookups_mongodb.get_transcript(db, transcript_id)
            print(transcript)
            print("rtanscript")
            variants_in_transcript = lookups_mongodb.get_variants_in_transcript(db, transcript_id)
            print("variants_in_transcript")
            coverage_stats = lookups_mongodb.get_coverage_for_bases(db, gene['xstart'] - EXON_PADDING, gene['xstop'] + EXON_PADDING)
            #add_transcript_coordinate_to_variants(db, variants_in_transcript, transcript_id)
            print("COVERAGE_STATS")
           # print(coverage_stats)
          
            # Add minicd info
            for i in range(0,len(variants_in_gene)):
                variants_in_gene[i]["minicd_info"] = "info_not_found"
                icd_code = variants_in_gene[i]["minicd"]
                icd10info = lookups_mongodb.get_icd_info(db, icd_code)

                if len(icd10info) > 0:


                    if "Name" in icd10info[0].keys():
                        variants_in_gene[i]["minicd_info"] = icd10info[0]["Name"]
            
            #print ("Transcripts in gene")
            #print(transcripts_in_gene)

            #print ("Variants in transcript")
            #print (variants_in_transcript)



            t = render_template(
                'gene_mongodb.html',
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


@app.route('/transcript/<transcript_id>')
def transcript_page(transcript_id):
    db = get_db()
    try:
        transcript = lookups_mongodb.get_transcript(db, transcript_id)

        cache_key = 't-transcript-{}'.format(transcript_id)
        t = cache.get(cache_key)
        if t is None:

            gene = lookups_mongodb.get_gene(db, transcript['gene_id'])
            gene['transcripts'] = lookups_mongodb.get_transcripts_in_gene(db, transcript['gene_id'])
            variants_in_transcript = lookups_mongodb.get_variants_in_transcript(db, transcript_id)

            coverage_stats = lookups_mongodb.get_coverage_for_transcript(db, transcript['xstart'] - EXON_PADDING, transcript['xstop'] + EXON_PADDING)

            add_transcript_coordinate_to_variants(db, variants_in_transcript, transcript_id)

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
            genes_in_region = lookups_mongodb.get_genes_in_region(db, chrom, start, stop)
            variants_in_region = lookups_mongodb.get_variants_in_region(db, chrom, start, stop)
            xstart = get_xpos(chrom, start)
            xstop = get_xpos(chrom, stop)
            coverage_array = lookups_mongodb.get_coverage_for_bases(db, xstart, xstop)
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
    db = get_db()
    try:
        variants = lookups_mongodb.get_variants_by_rsid(db, rsid)
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
        unsupported = "TTN" if query.upper() in lookups_mongodb.UNSUPPORTED_QUERIES else None
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
    datatype, identifier = lookups_mongodb.get_awesomebar_result(db, query)
    if datatype in ['gene', 'transcript']:
        gene = lookups_mongodb.get_gene(db, identifier)
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
    db = get_db()
    try:
        return render_template(
            'graph.html',
            graph_key=key
            )
    except Exception as e:
        print('Failed: %s' % e)
        abort(404)


@app.route('/meta/<key>')
def metapage(key):
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
    db = get_db()

    with open("./MRP_out/"+key+".mcmc.posteriors", "r") as inFile:
        probabilities = []
        admixture_data = []
        variants = []
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
            variants.append(variant)
           # print("variant",variant)
       # print("variants",variants)
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
           # print(len(gene_info))
            for j in range(1, int((len(gene_info)-1)/3)+1):
                gene["%s" % (j-1)] = gene_info[j]
            genes.append(gene)
        admixture_datagene.append(genes)
   # print("admixturedatagene",admixture_datagene)
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
#                headerarr.append(j)
                headeritem = header[j].decode('utf8').replace(" ","")
                headeritem = headeritem.rstrip(header[j][-3:]).upper()
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
                    headeritem = headeritem.rstrip(header[j][-3:]).upper()
                    tmp.append(headeritem.encode("ascii"))
                   # print(headeritem)
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
        return render_template(
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
    except Exception as e:
        print('Failed: %s' % e)
        abort(404)


@app.route('/login_page')
def login_page():
    return render_template('login_page.html')



#@basic_auth.required
@app.after_request
def apply_caching(response):
    # prevent click-jacking vulnerability identified by BITs
    response.headers["X-Frame-Options"] = "SAMEORIGIN"
    return response

if __name__ == "__main__":
    app.run(host = "0.0.0.0", port = 5000, debug=True, use_debugger=True, use_reloader=False)


