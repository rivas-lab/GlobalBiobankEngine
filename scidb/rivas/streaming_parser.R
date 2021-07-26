source('/home/scidb/GlobalBiobankEngine/scidb/database_load.R')

library(scidbstrm)
args = commandArgs(trailingOnly=TRUE)

if(args[1] == "parse_variants") 
{ 
  print("Calling parse variants")
  return(map(parse_variants))  
} else
{
  stop("Unknown arguments")
}