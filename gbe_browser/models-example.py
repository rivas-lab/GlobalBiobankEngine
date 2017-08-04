import numpy
import scidbpy

import lookups
import utils


db = scidbpy.connect()

#gene_names = ['IL33']          # Set to None to get all genes
gene_names = None
icds = ['HC49', 'HC382']        # Set to None to get all ICDs


# -- -
# -- - M-size Lists - --
# -- -

# SciDB lookup
# ---
gene_variant = lookups.get_gene_variant(db, gene_names=gene_names, icds=icds)

keys_M = []
for gv in gene_variant:
    # vep_annotations = lookups.parse_vep_annotations(gv['csq']['val'])
    # variant = {}
    # utils.add_consequence_to_variant(variant, vep_annotations)

    # major_consequence_M.append(variant['major_consequence'])
    # HGVSp_M.append(variant['HGVSp'])

    keys_M.append('{}-{}-{}-{}'.format(gv['chrom'],
                                       gv['pos'],
                                       gv['ref']['val'],
                                       gv['alt']['val']))

gene_names_M = gene_variant['gene_name']['val']
major_consequence_M = gene_variant['consequence']['val']
HGVSp_M = gene_variant['hgvsp']['val']


print(major_consequence_M)
print(HGVSp_M)
print(keys_M)
print(gene_names_M)


# -- -
# -- - M x N NumPy Arrays - --
# -- -

se_M_N = None
lor_M_N = None

for icd in icds:
    # SciDB lookup
    # ---
    variant_icd = lookups.get_variant_icd(db, gene_names=gene_names, icds=[icd])

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

print(lor_M_N)
print(se_M_N)
