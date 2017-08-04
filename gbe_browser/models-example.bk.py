import scidbpy

import lookups
import utils
import numpy

db = scidbpy.connect()

gene_names = ['RUNX3']          # Set to None to get all genes
icds = ['HC49', 'HC382']        # Set to None to get all ICDs


# SciDB lookups
# ---
# return 1-dimension NumPy arrays
gene_variant = lookups.get_gene_variant(db, gene_names=gene_names, icds=icds)
variant_icd = lookups.get_variant_icd(db, gene_names=gene_names, icds=icds)

print(gene_variant)
print(variant_icd)
print(len(gene_variant))
print(len(variant_icd))
#print(variant_icd[numpy.where(variant_icd['pos'] == 25243960)])
print(variant_icd[(variant_icd['icd']['val'] == "HC382") & (variant_icd['chrom'] == 1) & (variant_icd['pos'] == 25243960)])
# List of attributes available in each array
# ---
# Some attributes will have `null` and `val` sub-attributes. `null`
# stores the SciDB null code, while `val` stores the actual value (if
# attribute is not null)
print(gene_variant.dtype)
print(variant_icd.dtype)


# vep_annotations processing
# ---
# Extract vep_annotations of the first gene result
gene_id = gene_variant[0]['gene_id']['val']
csq = gene_variant[0]['csq']['val']
rsid = gene_variant[0]['rsid']
print("RSID")
print(rsid)
print('CHROM',gene_variant[0]['chrom'])
vep_annotations = lookups.parse_vep_annotations(csq, gene_id=gene_id)
variant = {}  # additional variant info
# Use exising utils function to populate extra variant info
utils.add_consequence_to_variant(variant, vep_annotations)

print(variant['category'])
print(variant['major_consequence'])
print(variant['HGVSp'])
