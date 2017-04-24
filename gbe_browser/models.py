"""
Greg McInes
Rivas Lab
gmcinnes@stanford.edu
January 15, 2017
"""

# Import packages
from __future__ import print_function
import argparse
import lookups
import numpy as np
import utils
import datetime
import subprocess
import os
from pymongo import MongoClient

"""
Main functional class
This class contains modeling functions used to generate genetic risk scores for ICD 10 codes from genotyping data.

Inputs:
    gene_names: A list of gene names to use
    icd_codes: a list of ICD 10 codes for which to build models
    threshold: p-value cutoff.  All betas with a p-value below the threshold will be returned
    lof: query loss of function variants only
    debug: print verbose output about operations
"""
class QueryGenome(object):
    def __init__(self, gene_names=None, icd_codes=None, threshold=1.0, category=None, all=False, debug=False):
        # Settings
        self.gene_names = gene_names
        self.icd_codes = icd_codes
        self.threshold = float(threshold)
        self.category = category
        self.all = all
        self.debug = debug

        # Set up database connection
        client = MongoClient('localhost', 27017)
        self.db = client.gbe

    """
    Execute the query genome function using preset values for gene_names and icd_codes
    """
    def run(self):
        return self.query_genome(self.gene_names, self.icd_codes)

    """
    Gather log odds ratios, standards errors, and annotation data for each variant as it relates to each ICD 10 code
    Input:
        gene_names: a list of gene names
        icd: a list of ICD 10 codes
    Return:
        Betas: A mxn numpy array containing log odd ratios for each variant and each ICD 10 code.  m=varaints, n=ICD codes
        Standard errors: A mxn numpy array containing standard errors for each variant and each ICD 10 code.
        Annotations: A list of variant annotations (variant consequence and category)
        Protein annotations: HGVSc values for each variant
        Key: chr-position-ref-alt
        Gene: A list of gene names associated with each variant
    """
    def query_genome(self, gene_names=None, icd=None):
        genes = []
        variant_count = 0
        icd_codes = []
        if not self.all:
            genes = list(self.get_genes(gene_names))
            for g in genes:
                self.get_gene_variants(g)
                if self.debug: print("%s variants added to gene %s" % (g.variant_count(), g.name))
                variant_count += g.variant_count()
                self.get_icd_data(g, icd)
                for v in g.variants:
                    icd_codes = utils.union(icd_codes, g.variants[v].icd.keys())
            print("%s total variant(s) found, %s ICD code(s)" % (variant_count, len(icd_codes)))
        else:
            # Get all variants, chromosome by chromosome.  Should reduce computation time.
            for c in range(1,23):
                xpos_start = int(str(c) + '000000000')
                xpos_end = int(str(c) + '999999999')
                chr = Gene('chr' + str(c), 'chr' + str(c), c, 0, 0, xpos_start, xpos_end)
                self.get_chromosome_variants(chr)
                if self.debug: print("%s variants added to chromosome %s" % (chr.variant_count(), chr.name))
                variant_count += chr.variant_count()
                self.get_icd_data(chr, icd)
                for v in chr.variants:
                    icd_codes = utils.union(icd_codes, chr.variants[v].icd.keys())
                genes.append(chr)
            print("%s total variant(s) found, %s ICD code(s)" % (variant_count, len(icd_codes)))
        return self.fill_matrices(genes, variant_count, icd_codes)

    """
    Populate return matrices for query_genome function
    Input:
        genes: A list of Gene objects
        variant_count: Total number of variants to be output, sum accross all genes
        icd_codes: A list of icd_codes to be included in output
    """
    def fill_matrices(self, genes, variant_count, icd_codes):
        icd_count = len(icd_codes)
        # Initialize arrays
        betas = np.zeros(shape=(variant_count, icd_count),dtype='float64')
        pvalues = np.zeros(shape=(variant_count, icd_count),dtype='float64')
        allele_frequencies = np.zeros(shape=(variant_count, 1),dtype='float64')
        se = np.zeros(shape=(variant_count, icd_count),dtype='float64')
        annotations = []
        protein_annotations = []
        variant_ids = []
        gene_return = []
        rsids = []
        icd = []
        alts = []
        icd_count = 0
        for i in icd_codes:
            variant_count = 0
            for g in genes:
                for v in g.variants:
                    if icd_count == 0:
                        variant_ids.append(str(g.variants[v].variant_id))
                        annotations.append(str(g.variants[v].annotations()))
                        protein_annotations.append(str(g.variants[v].hgvsc))
                        gene_return.append(g.name)
                        rsids.append(g.variants[v].rsid)
                        alts.append(g.variants[v].alt)
                        allele_frequencies[variant_count] = g.variants[v].allele_freq
                    if i in g.variants[v].icd and self.threshold > g.variants[v].icd[i].pvalue:
                        betas[variant_count][icd_count] = g.variants[v].icd[i].lor
                        se[variant_count][icd_count] = g.variants[v].icd[i].se
                        pvalues[variant_count][icd_count] = g.variants[v].icd[i].pvalue
                    else:
                        betas[variant_count][icd_count] = 0
                        se[variant_count][icd_count] = 0
                        pvalues[variant_count][icd_count] = 0
                    variant_count += 1
            icd.append(i)
            icd_count += 1
        #if self.debug:
        #    print("BETAS")
        #    print(betas)
        #    print("STANDARD ERROR")
        #    print(se)
        #    print("ANNOTATIONS")
        #    print(annotations)
        #    print("PROTEIN ANNOTATIONS")
        #    print(protein_annotations)
        #    print("VARIANT IDS")
        #    print(variant_ids)
        #    print("ICD")
        #    print(icd)
        #    print("GENES")
        #    print(gene_return)
        #    print("RSIDS")
        #    print(rsids)
        #    print("ALTS")
        #    print(alts)
        #    print("PVALUES")
        #    print(pvalues)
        #    print("ALLELE FREQUENCIES")
        #    print(allele_frequencies)
        return betas, se, pvalues, annotations, protein_annotations, \
               variant_ids, icd, gene_return, rsids, alts, allele_frequencies

    """
    A generator that yeilds Gene objects for each gene in gene_names. If gene_names is None, all genes will be returned
    Input:
        gene_names: A list of gene names
    Return:
        Yields Gene objects
    """
    def get_genes(self, gene_names):
        print(gene_names)
        if self.debug: print("Getting gene info")
        if not gene_names:
            if self.debug: print("No gene names passed.  Getting all genes.")
            result = self.db.genes.find(fields={'_id': False})
        else:
            if self.debug: print("Getting gene info for %s gene(s): %s" % (len(gene_names), ",".join(gene_names)))
            result = self.db.genes.find({'gene_name': {'$in': gene_names}}, fields={'_id': False}) # todo
        if result.count() == 0:
            print("No genes found!")
            exit()
        for r in result:
            new_gene = Gene(gene_id=r['gene_id'], name=r['gene_name'], chr=r['chrom'], start=r['start'], stop=r['stop'],
                            xstart=r['xstart'], xstop=r['xstop'])
            if self.debug: new_gene.print_gene()
            yield new_gene

    """
    Add Variant objects to Gene objects. This function queries the database for the given gene region and adds all found
    variants to the provided Gene object.
    Input:
        gene: A Gene object
    """
    def get_gene_variants(self, gene):
        raw_variants = lookups.get_variants_in_gene(self.db, gene.gene_id)
        for v in raw_variants:
            new_variant = Variant(v)
            if not self.category or new_variant.category in self.category:
                gene.variants[new_variant.xpos] = new_variant

    """
    Input:
        gene: A Gene object
    """
    def get_chromosome_variants(self, gene):
        raw_variants = lookups.get_variants_in_chromosome(self.db, str(gene.chr))
        for v in raw_variants:
            new_variant = Variant(v)
            if not self.category or new_variant.category in self.category:
                gene.variants[new_variant.xpos] = new_variant


    """
    Add all variants in the genome to a gene object
    """
    def get_all_variants(self, gene):
        raw_variants = lookups.get_all_variants(self.db)
        for v in raw_variants:
            new_variant = Variant(v)
            if not self.category or new_variant.category in self.category:
                gene.variants[new_variant.xpos] = new_variant

    """
    Add ICD objects to Variant objects, within a Gene object.  For each variant in a given Gene object, identify all
    ICD codes for which there is data in the database and add the log odds ratio and standard error to the associated
    Variant object.
    """
    def get_icd_data(self, gene, icd=None):
        if not icd:
            if self.debug: print("No ICD codes provided, fetching all.")
            results = self.db.icd.find({'xpos': {'$gte': gene.xstart, '$lte': gene.xstop}}, fields={'_id': False})
        else:
            if self.debug: print("Fetching ICD data for %s ICD code(s): %s" % (len(icd), ",".join(icd)))
            results = self.db.icd.find({'xpos': {'$gte': gene.xstart, '$lte': gene.xstop},
                                   'icd': {'$in': icd}},
                                  fields={'_id': False})
        for r in results:
            if r['xpos'] in gene.variants:
                new_icd = ICD(icd=r['icd'], pvalue=r['stats'][0]['pvalue'], lor=r['stats'][0]['lor'], se=r['stats'][0]['se'])
                gene.variants[r['xpos']].icd[r['icd']] = new_icd

"""
This Gene object stores all relevant information for each gene, including a dictionary of Variant objects.
Gene related functions should be added here.
"""
class Gene(object):
    def __init__(self, gene_id, name, chr, start, stop, xstart, xstop):
        self.gene_id = gene_id
        self.name = name
        self.chr = chr
        self.start = start
        self.stop = stop
        self.xstart = xstart
        self.xstop = xstop
        self.variants = {}

    """
    Print a brief summary of the gene
    """
    def print_gene(self):
        print(self.name, self.gene_id, self.chr, self.start, self.stop, len(self.variants))

    """
    Count the total number of variants stored within the gene
    """
    def variant_count(self):
        return len(self.variants)

"""
This Variant object stores all information relevant to each variant, including a dictionary of ICD objects.
The init function parses the raw response from the database query.
"""
class Variant(object):
    def __init__(self, variant):
        self.chr = None
        self.variant_id = None
        self.ref = None
        self.alt = None
        self.consequence = None
        self.category = None
        self.position = None
        self.rsid = None
        self.allele_freq = None
        self.xpos = None
        #self.vep_annotations = None
        self.hgvsc = []
        self.lof_info = []
        self.lof_filter = []
        self.parse_response(variant)
        self.icd = {}

    """
    Parse the response from the database and store each important value
    """
    def parse_response(self, v):
        self.chr = v['chrom']
        self.variant_id = v['variant_id']
        self.alt = v['alt']
        if 'major_consequence' in v:
            self.consequence = v['major_consequence']
        self.ref = v['ref']
        if 'category' in v:
            self.category = v['category']
        self.position = v['pos']
        self.rsid = v['rsid']
        self.allele_freq = v['allele_freq']
        self.xpos = v['xpos']
        if 'vep_annotations' in v:
            #self.vep_annotations = v['vep_annotations']
            for a in v['vep_annotations']:
                self.hgvsc.append(a['HGVSc'])
                if len(a['LoF_info']) > 0:
                    self.lof_info.append(a['LoF_info'])
                if len(a['LoF_filter']) > 0:
                    self.lof_filter.append(a['LoF_filter'])


    """
    Print a summary of the variant object
    """
    def print_variant(self):
        print("%s %s %s %s %s %s" % (self.variant_id, self.chr, self.position, self.ref, self.alt, self.consequence))

    """
    Return the total number of ICD codes stored for the variant
    """
    def icd_count(self):
        return len(self.icd)

    """
    Return the union of variant annotations
    """
    def annotations(self):
        return utils.union([self.consequence], [self.category])

"""
The ICD object stores information regarding how a variant relates to a specific ICD code
"""
class ICD(object):
    def __init__(self, icd, lor, se, pvalue):
        self.icd = icd
        self.lor = float(lor)
        self.se = float(se)
        self.pvalue = float(pvalue)

    def print_icd(self):
        print(self.icd, self.lor, self.se, self.pvalue)

# Questions for Manuel
# I'm using the alternate allele as the input allele for plink --score, is this the correct behavior?  Do you
# consider whether the reference or alternate was associated with the phenotype differently?
# I'm converting the nan's to 0's in the score file.  That should be fine, right?
# There are duplicate RSIDs returned with different betas
# The scores I'm getting all seem very low.  Would it be more valuable to only use a subset of variants?
# What about non-linear methods such as random forests?
# How come we don't use test and training sets?

# todo
# Create master score matrix to keep.  Save index of columns
class ComputeScore(object):
    def __init__(self, bfile=None, score_file=None, gene_names=None, icd_codes=None, output=None,
                 save_scores=False, threshold=1.0, category=None, all=False, snp_file=None, debug=False):
        self.gene_names = gene_names
        self.icd_codes = icd_codes
        self.bfile = bfile
        self.score_file = score_file
        self.file_name = output
        self.save_scores = save_scores
        self.threshold = float(threshold)
        self.category = category
        self.all = all
        self.snp_file = snp_file
        self.debug = debug
        self.to_cleanup = []
        self.snps_to_keep = []
        if not self.file_name:
            self.file_name = datetime.datetime.now().strftime("%Y%m%d%H%M%S")


    def run(self):
        # If the user has specified a file of snps to keep, read those into a dictionary
        if self.snp_file:
            self.process_snp_file()

        # Generate score file if it doesn't exist
        if not self.score_file:
            self.write_score_file()
            self.score_file = self.file_name + ".score"

        # Compute scores using plink
        if self.bfile:
            self.run_all_plink_scores()

            # Read the results and return them
            results = self.get_results()
            self.print_results(results)

        # Cleanup temp files
        self.cleanup()

    def process_snp_file(self):
        with open(self.snp_file) as f:
            for line in f:
                if line.startswith('#'): continue
                self.snps_to_keep.append(line.rstrip())

    def write_score_file(self):
        # Query genes from gene_names
        # Write rsid, alt, beta for each variant to file
        # Convert any nans to 0
        # If multiple icd codes, print a new column of scores for each icd code

        q = QueryGenome(gene_names=self.gene_names, icd_codes=self.icd_codes, threshold=self.threshold,
                        category=self.category, all=self.all, debug=self.debug)
        betas, se, pvalues, annotations, protein_annotations, variant_ids, icd, gene_return, rsids, alts, \
        allele_freqs = q.query_genome(gene_names=self.gene_names, icd=self.icd_codes)
        # Set self.icd to returned icd list to maintain order used in score file for later use
        self.icd_codes = icd
        if self.debug: print("Writing score file")

        f = open(self.file_name + ".score", 'w')
        # Create each line of output
        # RSID ALT ICD1SCORE ICD2SCORE etc.
        rsid_index = 0
        written_rsids = []
        for r in rsids:
            if r in written_rsids or r not in self.snps_to_keep:
                continue
            written_rsids.append(r)
            alt = alts[rsid_index]
            icd_index = 0
            outstring = "%s\t%s" % (r, alt)
            for i in icd:
                value = str(betas[rsid_index][icd_index])
                if value == "nan":
                    value = 0
                outstring += "\t%s" % value
                icd_index += 1
            print(outstring, file=f)
            rsid_index += 1
        f.close()
        if not self.save_scores:
            self.to_cleanup.append(self.file_name + ".score")
        else:
            f = open(self.file_name + ".icd", 'w')
            print("\n".join(icd), file=f)
            f.close()

    def run_all_plink_scores(self):
        # Need to work on this.  For this this assumes that the columns in the score file are in the same order
        # as the icd codes
        icd_index = 0
        for i in self.icd_codes:
            self.plink_score(i, icd_index + 3)
            icd_index += 1

    def plink_score(self, icd_code, score_column):
        if self.debug: print("Scoring %s" % icd_code)
        out = "%s-%s" % (self.file_name, icd_code)
        command = ['plink --score ' + self.score_file + ' 1 2 ' + str(score_column) +
                   ' --bfile ' + self.bfile + ' --out ' + str(out) + ' --noweb --allow-extra-chr']
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()
        if not process.returncode == 0:
            print("Plink command failed!\nPlease try again\n%s\n%s\n%s" % (" ".join(command), out, err))
        self.to_cleanup.append(out + ".profile")
        self.to_cleanup.append(out + ".log")
        self.to_cleanup.append(out + ".nopred")

    def get_results(self):
        results = {}
        for i in self.icd_codes:
            filename = "%s-%s.profile" % (self.file_name, i)
            line_count = 0
            with open(filename) as f:
                for line in f:
                    if line_count == 1:
                        fields = line.rstrip().split()
                        results[i] = fields[5]
                    line_count += 1
        return results

    def cleanup(self):
        for fn in self.to_cleanup:
            os.remove(fn) if os.path.exists(fn) else None

    def print_results(self, results):
        import operator
        sorted_by_value = sorted(results.items(), key=operator.itemgetter(1))
        for r in sorted_by_value:
            print("%s\t%s" % (r[0], r[1]))

class DiseaseDistributions(object):
    def __init__(self, icd, genes, debug=False):
        self.debug = debug
        self.icd = ['ICD1111']
        self.genes = genes

    def run(self):
        print("Running")
        # todo set p-value threshold, manny recommends 10-4
        q = QueryGenome(debug=True, threshold=0.0001)
        betas, se, pvalues, annotations, protein_annotations, \
        variant_ids, icd, gene_return, rsids, alts, allele_frequencies = q.query_genome(gene_names=self.genes, icd=self.icd)

        # Calculate mean
        # mean = log(p/1-p) + sum(2*f*log(r))
        # the second term is the sum of the allele frequencies multiplied by the log beta
        # Log beta value is what is stored in the plink .score file

        # Calculate variance
        # variance = sum(2*f*(1-f)*log(r)^2

        prevalence = np.array([0.0666680563256582])

        print(betas)
        print(allele_frequencies)
        print(prevalence)


        p_ratio = np.log10(np.divide(prevalence, np.subtract(1, prevalence)))  # todo, figure out if it is supposed to be log10
        risk_score = np.sum(np.multiply(2, np.multiply(betas, allele_frequencies)), axis=0)
        mean = p_ratio + risk_score

        variance = np.sum(np.multiply(2, np.multiply(allele_frequencies, np.multiply(np.subtract(1,allele_frequencies), np.multiply(betas, betas)))), axis=0)


        #variance = np.multiply((1np.multiply(betas, betas))

        print(p_ratio)
        print(mean)
        print(variance)





"""
Parse the command line
"""
def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'UK Biobank Modeling')
    parser.add_argument("-g", "--genes", nargs='+', type=str, help="Gene names to query.  Can include more than one.  "
                                                                   "If no gene names are provided, all genes will be "
                                                                   "queried. Ex: -g GENE_A GENE_B")
    parser.add_argument("-i", "--icd_codes", nargs='+', type=str, help="ICD 10 codes to query. If no ICD codes are "
                                                                       "provided, all ICD codes will be queried.")
    parser.add_argument("-b", "--bfile", help="Plink bfile for genetic risk scoring")
    parser.add_argument("--write_score_file", action='store_true', help="Save variant/ICD scores")
    parser.add_argument("--threshold", default=1.0, help="Value to threshold p-values of beta values, default 0")
    parser.add_argument("--category", nargs='+', type=str, help="Types of variants to allow")
    parser.add_argument("-d", "--debug", action='store_true', default=False,
                                help="Output debugging messages.  May be very verbose.")
    options = parser.parse_args()
    return options

"""
Main
"""
if __name__ == "__main__":
    options = parse_command_line()
    # todo check options here for which analysis to run
    #q = QueryGenome(gene_names=options.genes, icd_codes=options.icd_codes, debug=options.debug)
    #q.run()

    if options.write_score_file:
        c = ComputeScore(gene_names=options.genes, icd_codes=options.icd_codes, bfile=options.bfile,
                         save_scores=options.write_score_file, threshold=options.threshold, debug=options.debug)
        c.run()


    #d = DiseaseDistributions(genes=options.genes, icd=options.icd_codes)
    #d.run()
