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
library(MASS)
library(CompQuadForm)
DB = scidbconnect()

make_gv_query = function(gene_filter,
                         icd_filter,
                         variant_filter,
                         bim_filter)
{
  result = sprintf(
   "project(
      equi_join(
        apply(
          cross_join(
            project( icd, lor, se), 
            project( filter( icd_info, %s), Case),
            icd.icd_idx, icd_info.icd_idx
          ),
          icd_idx, icd_idx
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
          gene_idx, chrom, pos
        ),
        'left_names=chrom,pos',
        'right_names=chrom,pos',
        'algorithm=hash_replicate_right'
      ),
      gene_idx, icd_idx, chrom, pos, lor, se
    )", icd_filter, variant_filter, gene_filter, bim_filter)
  result
}

gene_var_lookup = function(gene_name      = "TP53", 
                           icds           = c("HC125","FH1065"),
                           variant_filter = "true", 
                           bim_filter     = "true")
{
  query = make_gv_query(gene_filter = sprintf("gene_name='%s'", gene_name),
                        icd_filter = paste(sprintf("icd='%s'", icds), collapse=" or "),
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


#Just sample functions to play with:
sample_function_mrp = function(gene_lookup_result)
{
library(MASS)
  #This returns an R named list - three numbers for each gene with text labels
  result = list()
#  write.table(gene_lookup_result, file = "test.1",quote = FALSE)
 # result$chrom = gene_lookup_result$chrom
 # result$pos = gene_lookup_result$pos
  result$chrom=gene_lookup_result[1,'chrom']
  result$pos=gene_lookup_result[1,'pos']
  result$nrows   = nrow(gene_lookup_result)
  b = gene_lookup_result$lor
  omegacor = matrix(c(1,.2284,-0.23911,-0.37297,0.2284,1,-0.08195,-0.01005,-0.23911,-0.08195,1,0.941319,-0.37297,-0.01005,0.941319,1),nrow = 4)
  thetacor = matrix(c(1,.159496,-.08871,-.14797,.159496,1,-.06311,-.08281,-.08871,-.06311,1,.873592,-.14797,-.08281,.873592,1),nrow = 4)
  
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
  d = eigen(nrows) - invA.sim %*% V
  d = d[d>1e-2]
  lbfsim = lbfsim + 0.5 * quad.form.sim
  result$l10bf = lbf/2.303
  result$l10bfsim = lbfsim/2.303
  result$pval = farebrother(quad.form.sim, d, maxit=100000, eps=10^-16, mode=1)$res
  return(result)
}

sample_function2 = function(gene_lookup_result)
{
  #This returns an R vector - of two - just two numbers for each gene
  a= nrow(gene_lookup_result)
  b= length(which(nchar(gene_lookup_result$ref) != 1 | nchar(gene_lookup_result$alt) != 1))
  return( c(a,b))
}

gene_var_stream = function(f, 
                           icds = c("HC125","FH1065"),
                           variant_filter = "true", 
                           bim_filter="true", 
                           genes)
{
  func_upload = as.scidb(DB, jsonlite::base64_enc(serialize(f, NULL)))
  gene_filter= 'true'
  if(!missing(genes))
  {
    gene_filter = paste0("gene_name='",genes,"'",collapse = " or ")
  }
  query = make_gv_query(gene_filter = gene_filter,
                        icd_filter = paste(sprintf("icd='%s'", icds), collapse=" or "),
                        variant_filter = variant_filter,
                        bim_filter = bim_filter)
  query = sprintf("
            faster_redimension(
              apply( %s, icd_idx_dbl, double(icd_idx), gene_idx_dbl, double(gene_idx), chrom_dbl, double(chrom), pos_dbl, double(pos)), 
              <gene_idx_dbl:double, icd_idx_dbl:double, chrom_dbl:double, pos_dbl:double, lor:double, se:double>
              [gene_idx = 0:*,1,0, n=0:*,10000000,0]
            )", query)
  
  query = sprintf("
   stream(
    %s,
    'Rscript /home/scidb/R_stream2/gene_icd_stream_child.R',
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
    project(filter(gene, %s), gene_name),
    'left_names=gene_idx', 
    'right_names=gene_idx',
    'right_outer=1'
   )", query, gene_filter)

  result =  iquery(DB, query, return=T)
  result =  data.frame(gene_name = result$gene_name, key=result$key, value=result$value, stringsAsFactors=F)
  result
}








