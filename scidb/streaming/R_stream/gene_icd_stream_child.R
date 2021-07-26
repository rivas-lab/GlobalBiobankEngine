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

library('scidbstrm')

function_to_apply = unserialize(jsonlite::base64_dec(getChunk()[[1]]))

process = function(input_frame)
{
  gene_idx = input_frame$gene_idx_dbl[1]  
  print(gene_idx)
  input = data.frame( input_frame[ , 2:ncol(input_frame) ])
  colnames(input) = c('chrom', 'pos', 'ref', 'alt', 'lor', 'se')

  success = T
  status = "OK"
  tryCatch({
    function_result = function_to_apply(input)
  }, error = function(e) {
    success <<- F
    status <<- paste("Error:", e)
  })
  if(success == F)
  {
    return(data.frame(gene_idx, key="status", value = status, stringsAsFactors=F))
  }
  if(class(function_result) == "logical" || class(function_result) == "character" || class(function_result) == "numeric" || 
     class(function_result) == "list" || class(function_result) == "integer")
  {
    if(is.null(names(function_result)))
    {
      key = rep("result", length(function_result))
    } else
    {
      key = names(function_result)
    }
    result = data.frame( 
      gene_idx   = gene_idx,
      key        = key,
      value      = as.character(function_result),
      stringsAsFactors=F
    )
    return(result)
  }
  result = data.frame(gene_idx = gene_idx, 
                      key="status", 
                      value = paste("Error: can't interpret result of class ", class(function_result)), stringsAsFactors=F)
  return(result)
}

print("running")

map(process)

