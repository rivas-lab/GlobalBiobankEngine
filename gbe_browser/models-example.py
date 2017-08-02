import scidbpy

import lookups
import utils


db = scidbpy.connect()

gene_names = ['RUNX3']          # Set to None to get all genes
icds = ['HC49', 'HC382']        # Set to None to get all ICDs

gene_names = ['SCYL3']
icds = ['INI6183']

## SciDB lookups
## ---
## return 1-dimension NumPy arrays
gene_variant = lookups.get_gene_variant(db, gene_names=gene_names, icds=icds)
variant_icd = lookups.get_variant_icd(db, gene_names=gene_names, icds=icds)

# print(len(gene_variant))
# print(len(variant_icd))

## List of attributes available in each array
## ---
## Some attributes will have `null` and `val` sub-attributes. `null`
## stores the SciDB null code, while `val` stores the actual value (if
## attribute is not null)
# print(gene_variant.dtype)
# print(variant_icd.dtype)


## vep_annotations processing
## ---
## Extract vep_annotations of the first gene result
# gene_id = gene_variant[0]['gene_id']['val']
# csq = gene_variant[0]['csq']['val']

# vep_annotations = lookups.parse_vep_annotations(csq, gene_id=gene_id)
# variant = {}  # additional variant info
## Use exising utils function to populate extra variant info
# utils.add_consequence_to_variant(variant, vep_annotations)

# print(variant['category'])
# print(variant['major_consequence'])
# print(variant['HGVSp'])

## Building size-M lists
## ---

major_consequence_M = []
HGVSp_M = []
keys_M = []
for gv in gene_variant:
    vep_annotations = lookups.parse_vep_annotations(gv['csq']['val'])
    variant = {}
    utils.add_consequence_to_variant(variant, vep_annotations)

    major_consequence_M.append(variant['major_consequence'])
    HGVSp_M.append(variant['HGVSp'])
    keys_M.append('{}-{}-{}-{}'.format(gv['chrom'],
                                       gv['pos'],
                                       gv['ref']['val'],
                                       gv['alt']['val']))
gene_names_M = gene_variant['gene_name']['val']

print(major_consequence_M)
print(HGVSp_M)
print(keys_M)
print(gene_names_M)
