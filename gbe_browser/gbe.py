from __future__ import division
import sys
import itertools
import json
import os
import pymongo
import pysam
import gzip
from parsing import *
import lookups
import random
import models
import numpy
import numpy.matlib
from targeted import mrpmm
from mr import mr
from graph import graph
from meta import meta
from utils import *
from flask import Flask, request, session, g, redirect, url_for, abort, render_template, flash, jsonify
from flask_compress import Compress
from flask import Flask, redirect, url_for, session
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

ADMINISTRATORS = (
    'gbe.help@gmail.com',
)

app = Flask(__name__)
app.debug=False

# Get Google API information
# client_secrets.json must be obtained from the Google Developers Console
with open('client_secrets.json') as secrets:
    secret_data = json.load(secrets)
google_client_id = secret_data['web']['client_id']
google_client_secret = secret_data['web']['client_secret']
redirect_uri = "/" + secret_data['web']['redirect_uris'][0].split('/')[-1]

SECRET_KEY = 'development key'
app.secret_key = SECRET_KEY
oauth = OAuth()

# Define google oauth protocol
google = oauth.remote_app('google',
                          base_url='https://www.google.com/accounts/',
                          authorize_url='https://accounts.google.com/o/oauth2/auth',
                          request_token_url=None,
                          request_token_params={'scope': 'https://www.googleapis.com/auth/userinfo.email',
                                                'response_type': 'code'},
                          access_token_url='https://accounts.google.com/o/oauth2/token',
                          access_token_method='POST',
                          access_token_params={'grant_type': 'authorization_code'},
                          consumer_key=google_client_id,
                          consumer_secret=google_client_secret)


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
    DB_HOST='mongodb',
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
    ICD_STATS_FILES=glob.glob(os.path.join(os.path.dirname(__file__), GBE_FILES_DIRECTORY, 'icdassoc','hybrid','c*.hybrid.rewrite.gz')),
    QT_STATS_FILES=glob.glob(os.path.join(os.path.dirname(__file__), GBE_FILES_DIRECTORY, 'icdassoc','hybrid','c*.linear.rewrite.gz')),
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
    client = pymongo.MongoClient(host=app.config['DB_HOST'], port=app.config['DB_PORT'])
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

def run_mrp(lof=True, missense=True, genes=None, fdr=5, phenidarr = ['ICD1462','ICD1463']):
    fdr = int(fdr)/100
    key = str(random.getrandbits(128))
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

        # Generate relevant files
        betas, se, pvalues, annotations, protein_annotations, variant_ids, icd, gene_return, rsids, alts, allele_frequencies = b.query_genome(genes,phenidarr)
        
        # Reshape betas from genome query
        C = numpy.matlib.eye(betas.shape[1])
        annotvec = [str(annotations[i].strip('"').strip('[').strip(']').strip("'")) for i in range(0,len(annotations))]

        # Run MRP with 2 clusters, output to the MRP_out subdirectory in the gbe_browser directory
        bicarr = []
        cmax = 6
        fail = 0
        for tmpc in range(1,cmax+1):
            [BIC, AIC, genedat] = mrpmm(betas,se, C, annotvec, gene_return, rsids, variant_ids,tmpc,key,C, numpy.linalg.inv(C),icd, fdr=fdr, niter=51,burn=10,thinning=1,verbose=True, outpath = './MRP_out/')
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
        if lbf > 1 and (lbf - lbf2) > 1 and lbf2 > 1:
#        if lbf > 1:
            [BIC, AIC, genedat] = mrpmm(betas,se, C, annotvec, gene_return, rsids, variant_ids,clustminval,key,C, numpy.linalg.inv(C), icd, fdr=fdr, niter=201,burn=100,thinning=1,verbose=True, outpath = './MRP_out/')
            print("log bayes factor: ",lbf)
            clustvalue = clustminval
        elif lbf2 > 1:
            clustminval = 2
            [BIC, AIC, genedat] = mrpmm(betas,se, C, annotvec, gene_return, rsids, variant_ids,clustminval,key,C, numpy.linalg.inv(C), icd, fdr=fdr, niter=201,burn=100,thinning=1,verbose=True, outpath = './MRP_out/')
            print("log bayes factor: ",lbf)
            clustvalue = clustminval
        else:
            [BIC, AIC, genedat] = mrpmm(betas,se, C, annotvec, gene_return, rsids, variant_ids,1,key,C, numpy.linalg.inv(C), icd, fdr=fdr, niter=51,burn=10,thinning=1,verbose=True, outpath = './MRP_out/')
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
    if check_credentials():
        return redirect(url_for('search_page'))
    return render_template('homepage.html')

@app.route('/search')
def search_page():
    if not check_credentials():
        return redirect(url_for('login'))
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
    if not check_credentials():
        return redirect(url_for('login'))
    db = get_db()
    try:
        chrom, pos = variant_str.split('-')
        pos = int(pos)
        # pos, ref, alt = get_minimal_representation(pos, ref, alt)
        xpos = get_xpos(chrom, pos)
        variant = lookups.get_variant(db, xpos)
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
        icdstats = lookups.get_variant_icd(db, xpos)
        indexes = []
        seend = {}
        for idx in range(0,len(icdstats)): 
            # ICD10=T81/Name=Complications of procedures, not elsewhere classified/Chapter=T/L95OR=0.97/U95OR=2.04/OR=1.40/pvalue=0.0756463/l10pval=1.12/Case=1229
            item = icdstats[idx]
            icd10 = item['icd']
            item['Code'] = icd10
            icd10info = lookups.get_icd_info(db, icd10)
            if len(icd10info) == 0:
                item['Name'] = 'NA'
                item['Group'] = 'NA'
                item['OR'] = 1
                item['L95OR'] = 1
                item['U95OR'] = 1
                item['pvalue'] = 1
                item['l10pval'] = 0
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
                elif icd10[0:3] == "MED":
                    item['Group'] = "MED"
                elif icd10[0:5] == "BRMRI":
                    item['Group'] = "BRMRI"
                else:
                    item['Group'] = icd10[0]
                item['OR'] = format(float(item['stats'][0]['or']), '.4g')
                item['L95OR'] = format(float(item['stats'][0]['l95or']), '.4g')
                item['U95OR'] = format(float(item['stats'][0]['u95or']), '.4g')
                item['pvalue'] = format(float(item['stats'][0]['pvalue']), '.4g')
                item['l10pval'] = format(float(item['stats'][0]['log10pvalue']), '.4g')
                item['Case'] = icd10info[0]['Case']
                se =  format(float(item['stats'][0]['se']), '.4g') 
                if float(item['l10pval']) <= 1 or float(se) >= .5 or int(item['Case']) <= 100 or item['Group'] == "INI" or item['Code'] == "HC67" or icd10 in seend:
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
            ordered_csqs=ordered_csqs
        )
    except Exception as e:
        print('Failed on variant:', variant_str, '; Error=', traceback.format_exc())
        abort(404)


@app.route('/coding/<icd_str>')
def icd_page(icd_str):
    if not check_credentials():
        return redirect(url_for('login'))
    db = get_db()
    try:
        icdlabel = icd_str
        #.strip('ADD').strip('INI').strip('BRMRI')
        # Try several cutoffs to see if there are too many variants to render.  Arbitrary 10k max.
        passing = False
        cutoff = None
        icd = None

        for p in [.001, .0001, .00001]:
            icd = lookups.get_icd_significant(db, str(icd_str), p)
           # print(icd_str,icd)
            if len(icd) < 100000:
                passing = True
            if passing:
                cutoff = p
               # print("CUTOFF",cutoff)
                break
        icd_info = lookups.get_icd_info(db, str(icdlabel))
        variants = get_icd_variant_table(icd_str, cutoff)
        print('Rendering ICD10: %s' % icdlabel)

        if len(icd_info) == 0:
            icd_info.append({'Case': 'NA', 'Name': 'NA', 'icd': icdlabel})
       # print(icd_info)


        return render_template(
            'icd.html',
            icd=icd,
            variants_in_gene = variants,
            icd_info=icd_info,
            cutoff=cutoff
            )
    except Exception as e:
        print('Failed on icd:', icd_str, '; Error=', traceback.format_exc())
        abort(404)


@app.route('/target/1')
def target_page():
    if not check_credentials():
        return redirect(url_for('login'))
    db = get_db()
    try:
        passing = False
        cutoff = None
        for p in [.00001]:
            vars = lookups.get_significant_prot(db, p, 0)
            print(vars)
            if len(vars) < 10000000000:
                passing = True
            if passing:
                cutoff = p
                break
        variants = get_prot_variant_table(0,cutoff)
        return render_template(
            'target.html',
            variants_in_gene = variants,
            cutoff=cutoff
            )
    except Exception as e:
        print('Failed on protective scan Error=', traceback.format_exc())
        abort(404)



def get_icd_variant_table(icd, p):
    db = get_db()
    significant_variant_ids = {}
    significant_variants = lookups.get_icd_significant(db, icd, p)
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
    if not check_credentials():
        return redirect(url_for('login'))
    db = get_db()
    try:
        print('Rendering Intensity page: %s' % affy_str)
        variant = lookups.get_variant_affy(db, affy_str)
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
                meta(key, betas, se, variant_ids,chains = 4, iter = 5000, warmup = 1000, cores = 2)
                return redirect('/meta/%s' % key)
        return get_gene_page_content(gene_id)




@app.route('/runMRP', methods=['GET', 'POST'])
def runMRP_page():

    if request.method == 'POST':

        genes = []

        missense = lof = False
        function = request.form['function']
        functionfdr = request.form['functionfdr']
        functionfdr = str(functionfdr)
        functionphen = request.form.getlist('phenotypes[]')
      #  print(functionphen)
        phenidarr = []
        for phenname in functionphen:
            phenidarr.append(str(phenname))
       # print(phenidarr)
        fdr = functionfdr.strip('[')
        fdr = fdr.strip(']')
       # print("fdrhtml",fdr)
       # print("functionfdr",functionfdr)
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
        #print(genes)
        # Run MRP function
        [key,lbf,clusterval, genedat] = run_mrp(lof=lof, missense=missense, genes=genes, fdr=fdr, phenidarr=phenidarr)
        # Send results of MRP function to mrp page

        return redirect('/mrp/%s' % key)


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
    if not check_credentials():
        return redirect(url_for('login'))
    db = get_db()
    try:
        gene = lookups.get_gene(db, gene_id)
        if gene is None:
            abort(404)
        cache_key = 't-gene-{}'.format(gene_id)
        t = cache.get(cache_key)
        if t is None:
           # print(gene_id)
            variants_in_gene = lookups.get_variants_in_gene(db, gene_id)
            transcripts_in_gene = lookups.get_transcripts_in_gene(db, gene_id)
            #print(variants_in_gene)
            # Get some canonical transcript and corresponding info
            transcript_id = gene['canonical_transcript']
            transcript = lookups.get_transcript(db, transcript_id)
            variants_in_transcript = lookups.get_variants_in_transcript(db, transcript_id)
            coverage_stats = lookups.get_coverage_for_transcript(db, transcript['xstart'] - EXON_PADDING, transcript['xstop'] + EXON_PADDING)
            add_transcript_coordinate_to_variants(db, variants_in_transcript, transcript_id)


            #print("\n\n\n\nVariants in gene")
            #print(variants_in_gene[0])

            # Add minicd info
            for i in range(0,len(variants_in_gene)):
                variants_in_gene[i]["minicd_info"] = "info_not_found"
                icd_code = variants_in_gene[i]["minicd"]
                icd10info = lookups.get_icd_info(db, icd_code)

                if len(icd10info) > 0:

                    
                    if "Name" in icd10info[0].keys():
                        names = icd10info[0]["Name"]
                        names = names.split()
                        names = "&nbsp;".join(names)
                        variants_in_gene[i]["minicd_info"] = names

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


@app.route('/transcript/<transcript_id>')
def transcript_page(transcript_id):
    if not check_credentials():
        return redirect(url_for('login'))
    db = get_db()
    try:
        transcript = lookups.get_transcript(db, transcript_id)

        cache_key = 't-transcript-{}'.format(transcript_id)
        t = cache.get(cache_key)
        if t is None:

            gene = lookups.get_gene(db, transcript['gene_id'])
            gene['transcripts'] = lookups.get_transcripts_in_gene(db, transcript['gene_id'])
            variants_in_transcript = lookups.get_variants_in_transcript(db, transcript_id)

            coverage_stats = lookups.get_coverage_for_transcript(db, transcript['xstart'] - EXON_PADDING, transcript['xstop'] + EXON_PADDING)

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
    if not check_credentials():
        return redirect(url_for('login'))
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
    if not check_credentials():
        return redirect(url_for('login'))
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
        gene = lookups.get_gene(db, identifier)
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

@app.route('/login')
def login():
    callback = url_for('authorized', _external=True)
    return google.authorize(callback=callback)

@app.route('/login_page')
def login_page():
    return render_template('login_page.html')


@app.route(redirect_uri)
@google.authorized_handler
def authorized(resp):
    access_token = resp['access_token']
    session['access_token'] = access_token, ''
    return redirect(url_for('search_page'))

@google.tokengetter
def get_acess_token():
    return session.get('access_token')

# Function to check Google Login Credentials for current session
def check_credentials():
    # Debugging.  Sick of logging in and getting redirected
    #if True:
    #    return True
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
        return True
    except URLError as e:
        if e.code == 401:
            # Unauthorized - bad token
            session.pop('access_token', None)
            return redirect(url_for('login_page'))
        return False


#@basic_auth.required
@app.after_request
def apply_caching(response):
    # prevent click-jacking vulnerability identified by BITs
    response.headers["X-Frame-Options"] = "SAMEORIGIN"
    return response

if __name__ == "__main__":
    app.run(host = "0.0.0.0", port = 5000, debug=False, use_debugger=False, use_reloader=False)

