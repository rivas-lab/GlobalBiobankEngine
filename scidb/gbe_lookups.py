import os
cwd = os.getcwd()

#first do this:
#cd ~
#git clone https://github.com/Paradigm4/SciDB-Py.git
#cd SciDB-Py
#git checkout devel 

os.chdir('/home/scidb/SciDB-Py')
from scidbpy.db import connect, iquery
from scidbpy.schema import Schema
os.chdir(cwd)

scidb = connect()
#It helps to get the schema of the ICD array, makes the query a little faster
icd_schema = iquery(scidb, "show(ICD)", fetch=True)
icd_schema=Schema.fromstring(icd_schema[0][0])

#Example lookup function for the ICDs. You specify chromosome, position, optionally the ICD tag:
# get_icd(1,795222, as_dataframe=True)
# get_icd(1,8999999,9999999, as_dataframe=True)
# get_icd(1,8999999,9999999,icd_string='1123', as_dataframe=True)
def get_icd(chromosome, start, end=None, icd_string=None, as_dataframe=False):
  if end is None:
    end = start
  query="between(ICD, null,"+ `chromosome`+","+`start`+",null,"\
                     "null,"+ `chromosome`+","+`end`+",null )"
  if icd_string is not None:
    query="project(cross_join("\
           "" + query + ","\
           "filter(ICD_INDEX, icd='"+icd_string+"'),"\
           "ICD.icd_id, ICD_INDEX.icd_id),"\
           "affyid, or_val, se, pvalue, lor, log10pvalue, l95or, u95or)"
  return iquery(scidb, query, fetch=True, as_dataframe=as_dataframe, schema=icd_schema)



