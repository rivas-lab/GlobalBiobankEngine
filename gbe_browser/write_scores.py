"""
Greg McInes
Rivas Lab
gmcinnes@stanford.edu
February 1, 2017
"""

import models
import argparse


def run(genes, icd_codes, threshold, category, all, snp_file, debug=False):
    c = models.ComputeScore(gene_names=genes, icd_codes=icd_codes,
                            save_scores=True, threshold=threshold,
                            category=category, all=all, snp_file=snp_file,
                            debug=debug)
    c.run()


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
    parser.add_argument("--threshold", default=1.0, help="Value to threshold p-values of beta values, default 1")
    parser.add_argument("--category", nargs='+', type=str, help="Types of variants to allow")
    parser.add_argument("--all", action='store_true', default=False, help="Fetch all variants available.  By default"
                                                                          "only coding variants will be retrieved.")
    parser.add_argument("--snp_file", default=None, help="Include only SNPs from provided file.  Must be new-line delimited rsids")
    parser.add_argument("-d", "--debug", action='store_true', default=False,
                                help="Output debugging messages.  May be very verbose.")
    options = parser.parse_args()
    return options

"""
Main
"""
if __name__ == "__main__":
    options = parse_command_line()
    run(genes=options.genes, icd_codes=options.icd_codes, threshold=options.threshold,
        category=options.category, all=options.all, snp_file=options.snp_file, debug=options.debug)

