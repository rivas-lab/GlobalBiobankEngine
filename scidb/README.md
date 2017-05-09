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
