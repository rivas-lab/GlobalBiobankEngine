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
.libPaths('/home/scidb/R/x86_64-pc-linux-gnu-library/3.4')
Sys.setenv("PKG_CXXFLAGS"="-I/usr/include")

function_to_apply = unserialize(jsonlite::base64_dec(getChunk()[[1]]))

process = function(input_frame)
{
  list_id = input_frame$list_id_dbl[1]
  icd_idx = input_frame$icd_idx_dbl[1]      
  colnames(input_frame) = c('list_id', 'icd_idx', 'icd', 'gene_name', 'chrom', 'pos', 'lor', 'se')
  success = T
  status = "OK"
  tryCatch({
    function_result = function_to_apply(input_frame)
  }, error = function(e) {
    success <<- F
    status <<- paste("Error:", e)
  })
  if(success == F)
  {
    return(data.frame(list_id,icd_idx, key="status", value = status, stringsAsFactors=F))
  }
  if(class(function_result) == "logical" || class(function_result) == "character" || 
     class(function_result) == "numeric" || 
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
      list_id   = list_id,
      icd_idx   = icd_idx,
      key        = key,
      value      = as.character(function_result),
      stringsAsFactors=F
    )
    print(result)
    return(result)
  }
  result = data.frame(list_id = list_id,
                      icd_idx = icd_idx, 
                      key="status", 
                      value = paste("Error: can't interpret result of class ", class(function_result)), stringsAsFactors=F)
  return(result)
}

map(process)

