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
import scidbpy

import config
import lookups

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
        self.db = scidbpy.connect()

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
        Protein annotations: HGVSp values for each variant
        Key: chr-position-ref-alt
        Gene: A list of gene names associated with each variant
    """
    def query_genome(self, gene_names=None, icd=None):
        genes = []
        icd_codes = []
        if not self.all:
            genes = gene_names
            gene_variant = lookups.get_gene_variant(self.db, gene_names = genes, icds = icd)
            results = lookups.get_variant_icd(self.db, gene_names=genes,icds=icd)
            icd_count = len(icd_codes)
            annotations = []
            protein_annotations = []
            variant_ids = []
            gene_return = []
            rsids = []
            icd = []
            alts = []
            icd_count = 0
            idxkeep = 0
            idxkeeparr = []
            idxo = open('idxo.txt','w')
            for v in gene_variant:
                chr = v['chrom']
                #        self.variant_id = v['variant_id']
                alt = v['alt']['val']
                if 'major_consequence' in v:
                    consequence = v['major_consequence']
                ref = v['ref']['val']
                if 'category' in v:
                    category = v['category']
                position = v['pos']
                rsid = v['rsid']['val']
                variant_id = str(chr) + "-" + str(position) + "-" + str(ref) + "-" + str(alt)
                allele_freq = v['exac_nfe']['val']  # !!!
                xpos = v['chrom'] * config.XOFF + v['pos']
                csq = v['csq']['val']
                gene_id = v['gene_id']['val']
                vep_annotations =  lookups.parse_vep_annotations(csq, gene_id=gene_id)
                variant = {}
                utils.add_consequence_to_variant(variant, vep_annotations)
                consequence = variant['major_consequence']
                category = variant['category']
                if category in self.category:
                    idxkeeparr.append(idxkeep)
                    idxo.write(str(idxkeep))
                idxkeep += 1
            gene_variant_subset = gene_variant[idxkeeparr]
            variant_countall = len(gene_variant_subset)
            # Initialize arrays
            betas = np.zeros(shape=(variant_countall, icd_count),dtype='float64')
            pvalues = np.zeros(shape=(variant_countall, icd_count),dtype='float64')
            allele_frequencies = np.zeros(shape=(variant_countall, 1),dtype='float64')
            se = np.zeros(shape=(variant_countall, icd_count),dtype='float64')
            variant_count = 0
            for v in gene_variant_subset:
                chr = v['chrom']
                #        self.variant_id = v['variant_id']
                alt = v['alt']['val']
                if 'major_consequence' in v:
                    consequence = v['major_consequence']
                ref = v['ref']['val']
                if 'category' in v:
                    category = v['category']
                position = v['pos']
                rsid = v['rsid']['val']
                variant_id = str(chr) + "-" + str(position) + "-" + str(ref) + "-" + str(alt)
                allele_freq = v['exac_nfe']['val']  # !!!
                idxo.write(variant_id)
                xpos = v['chrom'] * config.XOFF + v['pos']
                csq = v['csq']['val']
                gene_id = v['gene_id']['val']
                vep_annotations =  lookups.parse_vep_annotations(csq, gene_id=gene_id)
                variant = {}
                utils.add_consequence_to_variant(variant, vep_annotations)
                consequence = variant['major_consequence']
                category = variant['category']
                #        self.hgvsp.append(variant['HGVSp'])
                hgvsp = variant['HGVSp']
                if "LoF_info" in variant:
                    lof_info.append(variant['LoF_info'])
                if "LoF_filter" in variant:
                    lof_filter.append(variant['LoF_filter'])
                variant_ids.append(str(variant_id))
                annotations.append(str(consequence))
                protein_annotations.append(str(hgvsp))
                gene_return.append(gene_id)
                rsids.append(str(variant_id))
                alts.append(alt)
                allele_frequencies[variant_count,0] = allele_freq
                for i in icd_codes:
                    datareturn = variant_icd[(variant_icd['icd']['val'] == i) & (variant_icd['chrom'] == chr) & (variant_icd['pos'] == position)]
                    if len(datareturn) > 0:
                        betas[variant_count][icd_count] = datareturn['lor']['val']
                        se[variant_count][icd_count] = datareturn['se']['val']
                        pvalues[variant_count][icd_count] = datareturn['pvalue']['val']
                    else:
                        betas[variant_count][icd_count] = 0
                        se[variant_count][icd_count] = 0
                        pvalues[variant_count][icd_count] = 0
                    icd.append(i)
                    icd_count += 1
                variant_count += 1
            idxo.close()
        return betas, se, pvalues, annotations, protein_annotations, \
               variant_ids, icd, gene_return, rsids, alts, allele_frequencies


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
