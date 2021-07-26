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
.libPaths('/home/scidb/R/x86_64-pc-linux-gnu-library/3.4')
library(scidb)

DB = scidbconnect()

example_gene_list = c('EGFR', 'KRAS', 'PTEN', 'NOTCH2', 'NOTCH1', 'IL8')
example_gene_list2 = c('IL33', 'KRAS', 'IDH1', 'IDO1', 'PTEN', 'EGFR')
example_gene_lists = list(example_gene_list, example_gene_list2)

make_gene_list_frame = function(gene_lists)
{
  list_id = c()
  gene_name = c()
  n = 0
  for( i in 1:length(gene_lists))
  {
    sublist = gene_lists[[i]]
    list_id = c(list_id, rep(i, length(sublist)))
    gene_name = c(gene_name, sublist)    
  }
  return (data.frame(list_id, gene_name))
}

example_gene_list_frame = make_gene_list_frame(example_gene_lists)
gene_list_upload = NULL

make_gv_query = function(gene_list_frame,
                         icd_filter,
                         variant_filter,
                         bim_filter,
                         value_filter)
{
  gene_list_upload <<- as.scidb(DB, gene_list_frame, types=c('int64','string'))  
  gene_selection = sprintf("
    equi_join(
     project(
      apply(
       gene, gene_idx, gene_idx
      ),
      gene_idx, gene_name
     ),
     %s, 
     'left_names=gene_name',
     'right_names=gene_name'
    )", gene_list_upload@name)
  variant_selection = sprintf("
    project(
     equi_join(
      equi_join(
       project( filter( variant,%s), ref,alt),
       project(
        equi_join(
         variant_gene,
         %s,
         'left_names=gene_idx',
         'right_names=gene_idx',
         'algorithm=hash_replicate_right'
        ),
        list_id,
        gene_name,
        chrom,
        pos
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
     list_id, gene_name, chrom, pos, ref, alt
    )",
    variant_filter,
    gene_selection,
    bim_filter
   )  
  icd_selection = sprintf(
   "project(
      equi_join(
        apply(
          cross_join(
            project(
             filter(icd, %s),
             lor, se
            ), 
            project( filter( icd_info, %s), icd),
            icd.icd_idx, icd_info.icd_idx
          ),
          icd_idx, icd_idx
        ),
        project(
         %s,
         list_id, gene_name, chrom, pos
        ),
        'left_names=chrom,pos',
        'right_names=chrom,pos',
        'algorithm=hash_replicate_right'
      ),
      list_id, icd_idx, icd, gene_name, chrom, pos, lor, se
    )", value_filter, icd_filter, variant_selection)
  icd_selection
}

gene_var_lookup = function(gene_list     = c("TP53",'KRAS'),
                           icd           = c("HC125"),
                           variant_filter = "true", 
                           bim_filter     = "true",
                           value_filter   = "true")
{
  gene_lists = list(gene_list)
  gene_list_frame = make_gene_list_frame(gene_lists)
  query = make_gv_query(gene_list_frame,
                        icd_filter = sprintf("icd='%s'", icd),
                        variant_filter = variant_filter,
                        bim_filter = bim_filter,
                        value_filter = value_filter)
  result = iquery(DB, query, return=T)
  return(result[, 3:ncol(result)])
}

#Just sample functions to play with:
sample_function = function(gene_lookup_result)
{
  library(MASS)
 library(brms)
  #This returns an R named list - three numbers for each gene with text labels
  result = list()
  result$nrows   = nrow(gene_lookup_result)
  result$gene_name_str = paste(c(unique(gene_lookup_result$gene_name)),sep=",",collapse=",")
  result$n_genes = length(unique(gene_lookup_result$gene_name))
  result$icd     = gene_lookup_result$icd[1]
 result$chrom=gene_lookup_result[1,'chrom']
  result$pos=gene_lookup_result[1,'pos']
  b = gene_lookup_result$lor
  result$se = gene_lookup_result$se
  prior.V = diag(result$nrows)*.2*.2
  nrows = result$nrows
  Vsim = matrix(rep(1,nrows*nrows), nrows, nrows)*.9
  diag(Vsim) <- 1
  prior.Vsim = Vsim*.2*.2
  V = diag(result$nrows)*result$se*result$se
  A = prior.V+V
  A.sim = prior.Vsim + V
  Vdet = as.numeric(determinant(V,logarithm=TRUE)$modulus)
  lbf = -0.5*(as.numeric(determinant(A,logarithm=TRUE)$modulus)-Vdet)
  lbfsim = -0.5*(as.numeric(determinant(A.sim,logarithm=TRUE)$modulus)-Vdet)
  invA = ginv(A)
  invA.sim = ginv(A.sim)
  quad.form = t(b) %*% (ginv(V) - invA) %*% b
  lbf = lbf + 0.5 * quad.form
  quad.form.sim = t(b) %*% (ginv(V) - invA.sim) %*% b
  lbfsim = lbfsim + 0.5 * quad.form.sim
  result$l10bf = lbf/2.303
  result$l10bfsim = lbfsim/2.303
  result$genes = gene_lookup_result$gene_name 
  gsub("[\r\n]", "", gene_lookup_result$gene_name)
  gene_lookup_result$gene_name <- as.factor(gene_lookup_result$gene_name)
  levels(gene_lookup_result$gene_name) <- 1:length(levels(gene_lookup_result$gene_name))
  gene_lookup_result$gene_name <- as.numeric(gene_lookup_result$gene_name)
  result$gene_name = gene_lookup_result$gene_name
  posids =  as.factor(gene_lookup_result[,'pos'])
  levels(posids) <- 1:length(levels(posids))
  posids <- as.numeric(posids)
 # data <- list(b = b, ses = result$se, posids = posids, gene_name = gene_lookup_result$gene_name)
 # mod <- brm(b | se(ses, sigma = TRUE) ~ 1 + (1|posids) + gene_name ,  data = data, 
 #         cores = 4, control= list(adapt_delta=.99) )
 # result$summary <- summary(fit, waic = TRUE)
### Trying out this STAN stuff

  return(result)
}

gene_var_stream = function(f, 
                           icds = c("HC125","FH1065"),
                           variant_filter = "true", 
                           bim_filter="true", 
                           value_filter="true",
                           gene_lists = list(c('TP53', 'KRAS'),
                                             c('KRAS', 'PTEN', 'IDO1')))
{
  func_upload = as.scidb(DB, jsonlite::base64_enc(serialize(f, NULL)))
  gene_list_frame = make_gene_list_frame(gene_lists)
  icd_filter = paste(sprintf("icd='%s'", icds), collapse=" or ")
  query = make_gv_query(gene_list_frame,
                        icd_filter = icd_filter,
                        variant_filter = variant_filter,
                        bim_filter = bim_filter,
                        value_filter = value_filter)
  query = sprintf("
   faster_redimension(
    apply( %s, 
     icd_idx_dbl, double(icd_idx), 
     list_id_dbl, double(list_id), chrom_dbl, double(chrom), pos_dbl, double(pos)), 
     <list_id_dbl:double, icd_idx_dbl:double, icd:string, gene_name:string, chrom_dbl:double, 
      pos_dbl:double, lor:double, se:double>
     [list_id = 0:*,1,0, icd_idx = 0:*,1,0, n=0:*,10000000,0]
    )", query)
  
  query = sprintf("
   stream(
    %s,
    'Rscript /home/scidb/R_stream3/gene_icd_stream_child.R',
    'format=df',
    'types=double,double,string,string',
    'names=list_id_dbl,icd_idx_dbl,key,value',
    _sg(%s, 0)
   )",
   query,
   func_upload@name
  )
  query = sprintf("
   project(
    equi_join(
     apply(%s, list_id, int64(list_id_dbl),icd_idx, int64(icd_idx_dbl)),
     project(filter(icd_info, %s), icd),
     'left_names=icd_idx', 
     'right_names=icd_idx',
     'right_outer=1'
    ),
    list_id, icd_idx, icd, key, value 
   )", query, icd_filter)
  result =  iquery(DB, query, return=T, attributes_only=T)
  result = result[ ,3:ncol(result)]
  result = result[ order(result$list_id), ]
  result
}

