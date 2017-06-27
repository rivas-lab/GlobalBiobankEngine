# Loader

```bash
/home/scidb/GlobalBiobankEngine/scidb# python loader.py
make_fifo:FIFO:/tmp/tmpUrsUfL/fifo
make_fifo:FIFO:/tmp/tmpf1SNXq/fifo
make_fifo:FIFO:/tmp/tmpVxZ4hZ/fifo
make_fifo:FIFO:/tmp/tmp1J2Mu0/fifo
make_fifo:FIFO:/tmp/tmp0FfaXh/fifo
make_fifo:FIFO:/tmp/tmpPqf1my/fifo
make_fifo:FIFO:/tmp/tmpqFnd0C/fifo
make_fifo:FIFO:/tmp/tmp4FCpmY/fifo
load_qc:QC:file:/home/scidb/GlobalBiobankEngine/gbe_data/qc/UKBioBiLallfreqSNPexclude.dat lines:112 skip:1
load_qc:QC:file:/home/scidb/GlobalBiobankEngine/gbe_data/qc/ukb_ukbl_low_concordance.dat lines:417 skip:0
load_qc:Array:qc records:528
load_icd:Pipes:starting 8...
make_pipe:Spawn:zcat /home/scidb/GlobalBiobankEngine/gbe_data/icdassoc/hybrid/c1.RH56.phe.all.ff.PHENO1.glm.logistic.hybrid.rewrite.gz > /tmp/tmpUrsUfL/fifo pid:166
make_pipe:Spawn:zcat /home/scidb/GlobalBiobankEngine/gbe_data/icdassoc/hybrid/c1.RH116.phe.all.ff.PHENO1.glm.logistic.hybrid.rewrite.gz > /tmp/tmpf1SNXq/fifo pid:167
make_pipe:Spawn:zcat /home/scidb/GlobalBiobankEngine/gbe_data/icdassoc/hybrid/c1.RH29.phe.all.ff.PHENO1.glm.logistic.hybrid.rewrite.gz > /tmp/tmpVxZ4hZ/fifo pid:168
make_pipe:Spawn:zcat /home/scidb/GlobalBiobankEngine/gbe_data/icdassoc/hybrid/c1.HC410.phe.all.ff.PHENO1.glm.logistic.hybrid.rewrite.gz > /tmp/tmp1J2Mu0/fifo pid:169
make_pipe:Spawn:zcat /home/scidb/GlobalBiobankEngine/gbe_data/icdassoc/hybrid/c1.HC198.phe.all.ff.PHENO1.glm.logistic.hybrid.rewrite.gz > /tmp/tmp0FfaXh/fifo pid:170
make_pipe:Spawn:zcat /home/scidb/GlobalBiobankEngine/gbe_data/icdassoc/hybrid/c1.RH105.phe.all.ff.PHENO1.glm.logistic.hybrid.rewrite.gz > /tmp/tmpPqf1my/fifo pid:171
make_pipe:Spawn:zcat /home/scidb/GlobalBiobankEngine/gbe_data/icdassoc/hybrid/c1.HC366.phe.all.ff.PHENO1.glm.logistic.hybrid.rewrite.gz > /tmp/tmpqFnd0C/fifo pid:172
make_pipe:Spawn:zcat /home/scidb/GlobalBiobankEngine/gbe_data/icdassoc/hybrid/c1.RH7.phe.all.ff.PHENO1.glm.logistic.hybrid.rewrite.gz > /tmp/tmp4FCpmY/fifo pid:173
load_icd:Query:starting...
load_icd:Query:done
load_icd:Pipes:return code:0,0,0,0,0,0,0,0
load_icd:Pipes:starting 1...
make_pipe:Spawn:zcat /home/scidb/GlobalBiobankEngine/gbe_data/icdassoc/hybrid/c1.RH117.phe.all.ff.PHENO1.glm.logistic.hybrid.rewrite.gz > /tmp/tmpUrsUfL/fifo pid:182
load_icd:Query:starting...
load_icd:Query:done
load_icd:Pipes:return code:0
load_qt:Pipes:starting 2...
make_pipe:Spawn:zcat /home/scidb/GlobalBiobankEngine/gbe_data/icdassoc/hybrid/c1.6183.phe.coding.ff.all.initialdata.PHENO1.glm.linear.rewrite.gz > /tmp/tmpUrsUfL/fifo pid:184
make_pipe:Spawn:zcat /home/scidb/GlobalBiobankEngine/gbe_data/icdassoc/hybrid/c1.1389.phe.coding.ff.all.initialdata.PHENO1.glm.linear.rewrite.gz > /tmp/tmpf1SNXq/fifo pid:185
load_qt:Query:starting...
load_qt:Query:done
load_qt:Pipes:return code:0,0
remove_fifo:Remove:/tmp/tmpUrsUfL/fifo
remove_fifo:Remove:/tmp/tmpf1SNXq/fifo
remove_fifo:Remove:/tmp/tmpVxZ4hZ/fifo
remove_fifo:Remove:/tmp/tmp1J2Mu0/fifo
remove_fifo:Remove:/tmp/tmp0FfaXh/fifo
remove_fifo:Remove:/tmp/tmpPqf1my/fifo
remove_fifo:Remove:/tmp/tmpqFnd0C/fifo
remove_fifo:Remove:/tmp/tmp4FCpmY/fifo
```

# Lookups

```bash
/home/scidb/GlobalBiobankEngine/scidb# python lookups.py
/usr/local/lib/python2.7/dist-packages/scidbpy/schema.py:611: UserWarning: 9 type(s) promoted for null support. Precision loss may occur
  dtype=self.get_promo_atts_dtype())
   icd_id  chrom        pos  synthetic         affyid   or_val        se  \
0      24      1    4133241          0  Affx-10093612  3.74961  0.499301
1      24      1   39381448          0   Affx-9978429  7.12891  0.490708
2      24      1  208808979          0   Affx-6697267  3.94300  0.498864
3      24      2  125814441          0  Affx-36368492  8.50995  0.493669
4      24      2  129232604          0  Affx-17725509  4.71435  0.498569
5      24      2  153278499          0  Affx-18251591  4.03311  0.498311
6      24      2  195482149          0  Affx-18753501  4.06109  0.499543

     pvalue       lor  log10pvalue     l95or      u95or    icd
0  0.008121  1.321652     2.090399  1.409200   9.976994  RH117
1  0.000063  1.964158     4.203270  2.724734  18.651863  RH117
2  0.005957  1.371942     2.224963  1.483150  10.482585  RH117
3  0.000014  2.141236     4.841083  3.233758  22.394767  RH117
4  0.001870  1.550611     2.728133  1.774317  12.525998  RH117
5  0.005134  1.394538     2.289582  1.518690  10.710530  RH117
6  0.005024  1.401451     2.298915  1.525538  10.810909  RH117
```

---

# Loading the data

Currently the loader is in R. Source `database_load.R` and then do:
  1. recreate_db()
  2. load_dbsnp()
  3. load_variants()
  4. load_icd()
  5. load_icd_info()

There are some extraneous things like `make_extra_indeces` that aren't currently used, may speed things up in the future. Note the ICD loader currently loads the `hybrid` files only - which is the vast majority of the ICD data. Just to start with.

# Running ICD Lookups in Python

For a quick demo, first do this as user `scidb`:
```bash
$ cd ~
$ git clone https://github.com/Paradigm4/SciDB-Py.git
$ cd SciDB-Py
$ git checkout devel
```

Now you can do this from Python:
```python
>>> execfile('gbe_lookups.py')
>>> result = get_icd(1,795222) #return all ICDs for chromosome 1, position 795222
```

Note that `result` is a nested NumPy array - which is what the new SciDB-Py package returns by default. You can access the fields like this:
```python
>>> result[0]
((255, u'Affx-14445027'), (255, 0.599601), (255, 0.317647), (255, 0.107343), (255, -0.511490844976566), (255, 0.9692262713068632), (255, 0.321718298739391), (255, 1.1175036067570145), 0, 1, 795222, 0)
```

The dimensions, chromosome, position and synthetic are actually returned last at the moment:
```python
>>> result[0]["chrom"]
1
>>> result[0]["pos"]
795222
>>> result[0]["affyid"]
(255, u'Affx-14445027')
>>> result[0]["pvalue"]
(255, 0.107343)
```

Note that the fields `affyid`, `pvalue` (and all other attributes) are nullable by default. SciDB would allow them to be null. Thus they are returned as a tuple in the form of `(null code, value)` and null code 255 means "this value is not null". You could access the value portion like this:

```python
>>> result[0]["pvalue"][1]
0.10734299999999999
```

The Python package will also allow you to return the data as a Pandas Dataframe instead of a numpy array. This will automatically coerce nulls and present as a simpler structure. Use `as_dataframe=True` to do that:
```python
>>> result2 = get_icd(1,795222, as_dataframe=True)
>>> result2
           affyid    or_val        se    pvalue       lor  log10pvalue  \
0   Affx-14445027  0.599601  0.317647  0.107343 -0.511491     0.969226
1   Affx-14445027  0.971836  0.038178  0.454276 -0.028568     0.342680
2   Affx-14445027  1.010800  0.059992  0.857846  0.010742     0.066591
3   Affx-14445027  1.347090  0.262824  0.256948  0.297947     0.590155
4   Affx-14445027  1.192530  0.216579  0.416218  0.176077     0.380679
5   Affx-14445027  1.225790  0.135711  0.133572  0.203586     0.874285
6   Affx-14445027  0.870400  0.211242  0.511132 -0.138802     0.291467
7   Affx-14445027  1.058520  0.232931  0.807124  0.056872     0.093060
8   Affx-14445027  0.671942  0.216477  0.066268 -0.397583     1.178693
9   Affx-14445027  0.824208  0.405404  0.633443 -0.193332     0.198292
10  Affx-14445027  1.049170  0.265102  0.856321  0.047999     0.067363
11  Affx-14445027  1.137790  0.200774  0.520243  0.129088     0.283794
12  Affx-14445027  0.734139  0.320935  0.335553 -0.309057     0.474239

       l95or     u95or  icd_id  chrom     pos  synthetic
0   0.321718  1.117504       0      1  795222          0
1   0.901770  1.047347       1      1  795222          0
2   0.898668  1.136924       2      1  795222          0
3   0.804778  2.254846       3      1  795222          0
4   0.780034  1.823161       4      1  795222          0
5   0.939499  1.599321       5      1  795222          0
6   0.575315  1.316836       6      1  795222          0
7   0.670539  1.670990       7      1  795222          0
8   0.439605  1.027072       8      1  795222          0
9   0.372349  1.824415       9      1  795222          0
10  0.624003  1.764027      10      1  795222          0
11  0.767644  1.686415      11      1  795222          0
12  0.391375  1.377094      12      1  795222          0
```

Finally, there are a few options to try on the lookup function itself. For example:
```python
>>> get_icd(1,8999999,9999999) #look for a range instead of a specific coordinate
>>> get_icd(1,8999999,9999999,icd_string='1123') #restrict to a specific ICD name
```
