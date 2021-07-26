source('/home/scidb/GlobalBiobankEngine/scidb/database_load.R')
connect()

#Get a sense of what data is loaded
as.R(DB$project(DB$list(), name)) #This is the more user-friendly form DB$ will auto-complete. DB is set in connect()
iquery(DB, "project(list(), name)", return=T) #Same thing

as.R(DB$summarize(ICD)) #stats: count size and chunks
as.R(DB$limit(ICD, 10)) #a selection of 10 entries
as.R(DB$summarize(VARIANT)) #do it with another array
as.R(DB$limit(VARIANT, 10)) #a selection of 10 entries

#ICD lookup
#get_icd(1, 795222);                                        #position
#get_icd(1, 795222, icd_string='1126');                     #position, pheno
#get_icd(1, start=795222, end=1795222);                     #range
#get_icd(1, start=795222, end=1795222, icd_string='1126');  #range, pheno
get_icd = function(chrom, start, end, icd_string)
{
  if(missing(end))
  {
    end=start
  }
  if(missing(icd_string))
  {
    iquery(DB,
     sprintf(
      "between(ICD, null, %i, %i, null,  null, %i, %i, null)",
      chrom,start, chrom, end 
     ),
     return=T,
     schema=ARRAYS[["ICD"]]
    )
  }
  else
  {
    iquery(DB,
     sprintf(
       "project(
         cross_join(
          between(ICD, null, %i, %i, null,  null, %i, %i, null),
          filter(ICD_INDEX, icd='%s'),
          ICD.icd_id, ICD_INDEX.icd_id
         ),
         affyid, or_val, se, pvalue, lor, log10pvalue, l95or, u95or
       )",
       chrom,start, chrom, end, icd_string
     ),
     return=T,
     schema=ARRAYS[["ICD"]]
    )
  }
}

#Variant lookup
#get_variant(1,889238) 
#x= get_variant(1,889238, 890238) #careful - printing a lot of these can cause trouble for rstudio
get_variant = function(chrom, start, end)
{
  if(missing(end))
  {
    end=start
  }
  iquery(DB,
   sprintf(
    "between(VARIANT, %i, %i, null, %i, %i, null)",
    chrom,start, chrom, end 
   ),
   return=T,
   schema=ARRAYS[["VARIANT"]]
  )
}

## ## Some prototypes of Manny's query below:
challenge_query = function()
{
  t1=proc.time();
  a = iquery(DB,
    "project(
      apply(
       cross_join(
        cross_join(
         project(filter(ICD, (lor<0 or or_val<1) and pvalue < 0.001), lor, or_val, pvalue),
         project(filter(ICD_INFO, case>300), icd_name, case),
         ICD.icd_id,
         ICD_INFO.icd_id
        ),
        project(filter(VARIANT, regex(info, '.*stop_gained.*')), ref, alt),
        ICD.chrom, 
        VARIANT.chrom,
        ICD.pos, 
        VARIANT.pos
       ),
      chromosome, chrom, position, pos
     ),
     chromosome, position, ref, alt, lor, or_val, pvalue, icd_name, case
    )",
     return=T,
     only_attributes=T
  )
  proc.time()-t1
  head(a)
  nrow(a)
}

challenge_query2 = function()
{
  t1=proc.time();
  a = iquery(DB,
    "project(
      apply(
       cross_join(
        cross_join(
         project(filter(ICD, (lor<0 or or_val<1) and pvalue < 0.001), lor, or_val, pvalue),
         project(filter(ICD_INFO, case>300), icd_name, case),
         ICD.icd_id,
         ICD_INFO.icd_id
        ),
        project(
         filter(
          apply(VARIANT, effect, iif(regex(info, '.*stop_gained.*'),       'stop_gained',
                                 iif(regex(info, '.*splice_acceptor.*'),   'splice_acceptor',
                                 iif(regex(info, '.*splice_donor.*'),      'splice_donor',
                                 iif(regex(info, '.*frameshift.*'),        'frameshift', 
                                                                           null))))
          ),
          effect is not null
         ), 
         ref, alt, effect
        ),
        ICD.chrom, 
        VARIANT.chrom,
        ICD.pos, 
        VARIANT.pos
       ),
      chromosome, chrom, position, pos
     ),
     chromosome, position, ref, alt, effect, lor, or_val, pvalue, icd_name, case
    )",
             return=T,
             only_attributes=T
  )
  proc.time()-t1
  head(a)
  nrow(a)
}