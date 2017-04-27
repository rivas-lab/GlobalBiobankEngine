source('/home/scidb/GlobalBiobankEngine/scidb/database_load.R')
connect()






#get_icd(1, 795222);                                        0.1s
#get_icd(1, 795222, icd_string='1387');                     0.2s
#get_icd(1, start=795222, end=1795222);                     0.37s
#get_icd(1, start=795222, end=1795222, icd_string='1387');  0.37s
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
         affyid, or_val, se, pvalue
       )",
       chrom,start, chrom, end, icd_string
     ),
     return=T,
     schema=ARRAYS[["ICD"]]
    )
  }
}

#get_variant(1,889238) : 0.12s
#x= get_variant(1,889238, 1889238) : 0.16s
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