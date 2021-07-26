#BEGIN_COPYRIGHT
#
# PARADIGM4 INC.
# This file is distributed along with the Paradigm4 Enterprise SciDB 
# distribution kit and may only be used with a valid Paradigm4 contract 
# and in accord with the terms and conditions specified by that contract.
#
# Copyright (C) 2010 - 2016 Paradigm4 Inc.
# All Rights Reserved.
#
#END_COPYRIGHT

library(scidb)

DB = scidbconnect()

make_gv_query = function(gene_filter,
                         icd,
                         variant_filter,
                         bim_filter)
{
  result = sprintf(
   "project(
      equi_join(
        cross_join(
          project( icd, lor, se), 
          project( filter( icd_info, icd='%s'), Case),
          icd.icd_idx, icd_info.icd_idx
        ),
        project(
          equi_join(
            equi_join(
              project( filter( variant,%s), rsid, ref, alt, exac_nfe),
              project(
                equi_join(
                  variant_gene,
                  equi_join(
                    gene_index,
                    project( filter( gene, %s), gene_name),
                    'left_names=gene_idx',
                    'right_names=gene_idx',
                    'algorithm=hash_replicate_right'
                  ),
                  'left_names=gene_idx',
                  'right_names=gene_idx',
                  'keep_dimensions=1',
                  'algorithm=hash_replicate_right'
                ),
                chrom,
                pos,
                gene_idx
              ),
              'left_names=chrom,pos',
              'right_names=chrom,pos',
              'algorithm=hash_replicate_right' 
            ),
            project( filter(bim, %s), ref, alt ),
            'left_names=chrom,pos,ref,alt',
            'right_names=chrom,pos,ref,alt',
            'algorithm=hash_replicate_left'
          ),
          gene_idx, chrom, pos, ref, alt
        ),
        'left_names=chrom,pos',
        'right_names=chrom,pos',
        'algorithm=hash_replicate_right'
      ),
      gene_idx, chrom, pos, ref, alt, lor, se
    )", icd, variant_filter, gene_filter, bim_filter)
  result
}

gene_var_lookup = function(gene_name      = "TP53", 
                           icd            = "HC125",
                           variant_filter = "true", 
                           bim_filter     = "true")
{
  query = make_gv_query(gene_filter = sprintf("gene_name='%s'", gene_name),
                        icd = icd,
                        variant_filter = variant_filter,
                        bim_filter = bim_filter)
  result = iquery(DB, query, return=T)
  return(result[, 4:ncol(result)])
}

#Just sample functions to play with:
sample_function = function(gene_lookup_result)
{
  #This returns an R named list - three numbers for each gene with text labels
  result = list()
  result$nrows   = nrow(gene_lookup_result)
  result$min_lor = min(gene_lookup_result$lor)
  result$avg_se = mean(gene_lookup_result$se)
  return(result)
}

sample_function2 = function(gene_lookup_result)
{
  #This returns an R vector - of two - just two numbers for each gene
  a= nrow(gene_lookup_result)
  b= length(which(nchar(gene_lookup_result$ref) != 1 | nchar(gene_lookup_result$alt) != 1))
  return( c(a,b))
}

gene_var_stream = function(f, icd = "HC125", variant_filter = "true", bim_filter="true", genes)
{
  func_upload = as.scidb(DB, jsonlite::base64_enc(serialize(f, NULL)))
  gene_filter= 'true'
  if(!missing(genes))
  {
    gene_filter = paste0("gene_name='",genes,"'",collapse = " or ")
  }
  query = make_gv_query(gene_filter = gene_filter,
                        icd = icd,
                        variant_filter = variant_filter,
                        bim_filter = bim_filter)
  query = sprintf("
            faster_redimension(
              apply( %s, gene_idx_dbl, double(gene_idx), chrom_dbl, double(chrom), pos_dbl, double(pos)), 
              <gene_idx_dbl:double, chrom_dbl:double, pos_dbl:double, ref:string, alt:string, lor:double, se:double>
              [gene_idx = 0:*,1,0, n=0:*,10000000,0]
            )", query)
  
  query = sprintf("
   stream(
    %s,
    'Rscript /home/scidb/R_stream/gene_icd_stream_child.R',
    'format=df',
    'types=double,string,string',
    'names=gene_idx_dbl,key,value',
    _sg(%s, 0)
   )",
   query,
   func_upload@name
  )
  query = sprintf("
   equi_join(
    apply(%s, gene_idx, int64(gene_idx_dbl)),
    project(gene, gene_name),
    'left_names=gene_idx', 
    'right_names=gene_idx',
    'right_outer=1'
   )", query)

  result =  iquery(DB, query, return=T)
  result =  data.frame(gene_name = result$gene_name, key=result$key, value=result$value, stringsAsFactors=F)
  result
}








