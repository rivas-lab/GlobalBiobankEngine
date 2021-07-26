
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
                         icds,
                         variant_filter,
                         icd_filter,
                         bim_filter)
{
  result = sprintf(
   "project(
      equi_join(
        apply(
          cross_join(
            project( filter(icd, %s), lor, se), 
            project( filter( icd_info, %s), icd),
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
            project( filter(bim, %s), ref, alt, consequence ),
            'left_names=chrom,pos,ref,alt',
            'right_names=chrom,pos,ref,alt',
            'algorithm=hash_replicate_left'
          ),
          gene_idx, chrom, pos, consequence
        ),
        'left_names=chrom,pos',
        'right_names=chrom,pos',
        'algorithm=hash_replicate_right'
      ),
      gene_idx, icd_idx, icd, chrom, pos, lor, se, consequence
    )",  icd_filter, icds, variant_filter, gene_filter, bim_filter)
  result
}

gene_var_lookup = function(gene_name      = "TP53", 
                           icds           = c("HC125","FH1065"),
                           variant_filter = "true", 
                           bim_filter     = "true",
                           icd_filter     = "true")
{
  query = make_gv_query(gene_filter = sprintf("gene_name='%s'", gene_name),
                        icds = paste(sprintf("icd='%s'", icds), collapse=" or "),
                        variant_filter = variant_filter,
                        bim_filter = bim_filter,
                        icd_filter = icd_filter)
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
  result$first_consequence = gene_lookup_result$consequence[1]
  return(result)
}


#Just sample functions to play with:
sample_function_mrp = function(gene_lookup_result)
{
library(MASS)
library(CompQuadForm)
  #This returns an R named list - three numbers for each gene with text labels
  result = list()
#  write.table(gene_lookup_result, file = "test.1",quote = FALSE)
 # result$chrom = gene_lookup_result$chrom
 # result$pos = gene_lookup_result$pos
  result$chrom=gene_lookup_result[1,'chrom']
  result$pos=gene_lookup_result[1,'pos']
  result$icd_1=min(gene_lookup_result$icd)
  result$icd_2=max(gene_lookup_result$icd)
  result$nrows   = nrow(gene_lookup_result)
  b = gene_lookup_result$lor
  bprot = b
  result$csq = gene_lookup_result$consequence
   icds <- list()
# glaucoma
#  icds['HC276'] = 1 
#  icds['INI5263'] = 2 
  icds['INI5255'] = 1 
# asthma
#  icds['HC382'] = 1 
#  icds['INI30150'] = 2 
#  icds['INI3062'] = 3 
#  icds['INI3063'] = 4
  omegacor = matrix(c(1,.2284,-0.23911,-0.37297,0.2284,1,-0.08195,-0.01005,-0.23911,-0.08195,1,0.941319,-0.37297,-0.01005,0.941319,1),nrow = 4)
  omegacor = matrix(c(1,.66,.27,.66,1,.80, .27, .80, 1),nrow = 3)
  thetacor = matrix(c(1,.05,.013,.05,1,.57, .013, .57, 1),nrow = 3)
#  omegacor <- c(1)
#  thetacor <- as.matrix(c(1))
#  thetacor = matrix(c(1,.159496,-.08871,-.14797,.159496,1,-.06311,-.08281,-.08871,-.06311,1,.873592,-.14797,-.08281,.873592,1),nrow = 4)
# thetacor = omegacor
  result$se = gene_lookup_result$se
  sd.lof = .5
  sd.mis = .2
  prior.V = diag(result$nrows)*sd.mis*sd.mis
  result$nr = dim(b)[1]
  result$nc = length(b)
  nrows = result$nrows
  Vsim = matrix(rep(1,nrows*nrows), nrows, nrows)*.99999
  diag(Vsim) <- 1
  prior.Vsim = Vsim*sd.mis*sd.mis
  lof.annot = c("stop_gained","frameshift_variant","splice_acceptor_variant","splice_donor_variant")
  V = diag(result$nrows)*result$se*result$se
  # result$test = icds[[gene_lookup_result[,'icd']]]
  result$test2 = gene_lookup_result$icd
  for(i in 1:result$nrows){
     if(icds[[gene_lookup_result[i,'icd']]] < 2){
        if(result$csq[i] %in% lof.annot){
      #   bprot[i] = bprot[i] + .5
     }
      else{
       #     bprot[i] = bprot[i] + .2
          }
      }
else{
            bprot[i] = bprot[i] 
  }
      for(j in 1:result$nrows){
          if(gene_lookup_result[i,'pos'] == gene_lookup_result[j,'pos']){
                  V[i,j] = result$se[i]*result$se[j]*thetacor[icds[[gene_lookup_result[i,'icd']]],icds[[gene_lookup_result[j,'icd']]]]
               prior.V[i,j] = sd.mis*sd.mis*omegacor[icds[[gene_lookup_result[i,'icd']]],icds[[gene_lookup_result[j,'icd']]]]
               }
                  if((result$csq[i] %in% lof.annot) && !(result$csq[j] %in% lof.annot) ){
                              prior.Vsim[i,j] = sd.mis*sd.lof*omegacor[icds[[gene_lookup_result[i,'icd']]],icds[[gene_lookup_result[j,'icd']]]]
       
                      }
           else if((result$csq[j] %in% lof.annot) && !(result$csq[i] %in% lof.annot) ){
	   	             prior.Vsim[i,j] = sd.mis*sd.lof*omegacor[icds[[gene_lookup_result[i,'icd']]],icds[[gene_lookup_result[j,'icd']]]]
       		      }
           else if((result$csq[i] %in% lof.annot) && (result$csq[j] %in% lof.annot) ){
	   	             prior.Vsim[i,j] = sd.lof*sd.lof*omegacor[icds[[gene_lookup_result[i,'icd']]],icds[[gene_lookup_result[j,'icd']]]]
          if(gene_lookup_result[i,'pos'] == gene_lookup_result[j,'pos']){                       
      prior.V[i,j] = sd.lof*sd.lof*omegacor[icds[[gene_lookup_result[i,'icd']]],icds[[gene_lookup_result[j,'icd']]]]
                      }
		      }
           }
    }
                
        
  A = prior.V+V
  A.sim = prior.Vsim + V
  Vdet = as.numeric(determinant(V,logarithm=TRUE)$modulus)
  lbf = -0.5*(as.numeric(determinant(A,logarithm=TRUE)$modulus)-Vdet)
  lbfsim = -0.5*(as.numeric(determinant(A.sim,logarithm=TRUE)$modulus)-Vdet)
  invA = ginv(A)
  invA.sim = ginv(A.sim)
  quad.form = t(b) %*% (ginv(V) - invA) %*% b
# quad.form = t(b) %*% ginv(V) %*% b - t(bprot) %*% invA %*% bprot
  lbf = lbf + 0.5 * quad.form
  quad.form.sim = t(b) %*% (ginv(V) - invA.sim) %*% b
# quad.form.sim = t(b) %*% ginv(V) %*% b - t(bprot) %*% invA.sim %*% bprot
#  lbfsim = lbfsim + 0.5 * quad.form.sim
#  result$l10bf = lbf/2.303
#  result$l10bfsim = lbfsim/2.303
   B = diag(nrows) - invA.sim %*% V
  d = eigen(B)$values
  d = d[d>1e-2]
  lbfsim = lbfsim + 0.5 * quad.form.sim
  result$l10bf = lbf/2.303
  result$l10bfsim = lbfsim/2.303
  result$pvalue = 0 
  result$pvalue = farebrother(quad.form.sim, d, maxit=100000, eps=10^-16, mode=1)$Qq

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
                           icd_filter="true", 
                           genes)
{
  func_upload = as.scidb(DB, jsonlite::base64_enc(serialize(f, NULL)))
  gene_filter= 'true'
  if(!missing(genes))
  {
    gene_filter = paste0("gene_name='",genes,"'",collapse = " or ")
  }
  query = make_gv_query(gene_filter = gene_filter,
                        icds = paste(sprintf("icd='%s'", icds), collapse=" or "),
                        variant_filter = variant_filter,
                        bim_filter = bim_filter, 
                        icd_filter = icd_filter)
  query = sprintf("
            faster_redimension(
              apply( %s, icd_idx_dbl, double(icd_idx), gene_idx_dbl, double(gene_idx), chrom_dbl, double(chrom), pos_dbl, double(pos)), 
              <gene_idx_dbl:double, icd_idx_dbl:double, icd:string, chrom_dbl:double, pos_dbl:double, lor:double, se:double, consequence:string>
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


