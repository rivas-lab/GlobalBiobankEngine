# Formatting note for icdinfo:

`icdinfo.txt` contains six columns for each entry, which should be formatted as below. An example follows:
```icd_code	n_cases	phenotype_description	n_cases	n_cases	Y
BIN20483	2433	Ever_attempted_suicide	2433	2433	Y
```
Note that columns 2, 4, and 5 contain duplicate information, and that all lines will contain the terminating "Y" column irrespective of their contents. Summary statistic files uploaded to either GBE server (production or dev) should be named according to some standard (which may vary by project), but that contains the icd code in the filename. For brevity and ease of encoding, the icd code should *not* be the same as the descriptor.
