library(scidb)
DB = NULL

source('/home/scidb/GlobalBiobankEngine/scidb/constants.R')
ARRAYS=list()

connect = function()
{
  if(is.null(DB))
  {
    DB <<- scidbconnect()
  }
}

recreate_db = function()
{
  confirm = readline(prompt = "This will erase the DB, are you sure? (y/n): ")
  connect()
  if(confirm != 'y' && confirm != 'Y')
  {
    print("Exiting")
    return();
  }
  for( i in 1:length(ARRAYS))
  {
    name = names(ARRAYS)[i]
    schema = ARRAYS[[i]]
    tryCatch({ iquery(DB, paste0("remove(", name, ")"))  },  warning = invisible, error = invisible )
    iquery(DB, paste0("create array ", name, " ", schema))
  }
}

remove_versions = function(array)
{
  connect()
  v = iquery(DB, sprintf("versions(%s)", array), return=TRUE)
  if(nrow(v) == 0)
  {
    return();
  }
  mv = max(v$version_id)
  iquery(DB, sprintf("remove_versions(%s, %i)", array, mv))
}

run_pipe = function(file, command = 'zcat', path = '/tmp/fifo')
{
  system(paste("rm -rf", path))
  system(paste("mkfifo", path))  
  system(paste(command, file, ">", path), wait=FALSE)
}

ARRAYS[["DBSNP"]] = "<rsid:int64> [chrom=1:25,1,0, pos = 0:*,10000000,0, synthetic=0:999,1000,0]"

load_dbsnp = function()
{
  connect()
  run_pipe(DBSNP_FILE, path='/tmp/fifo', command='zcat')
  iquery(DB, "
   insert(
     redimension(
       apply(
        filter(
          aio_input('/tmp/fifo', 'num_attributes=4'),
          a1 <> 'PAR'
        ),
        rsid,  int64(a0),
        chrom, int64(rsub(a1, 's/T+$//g')),
        pos,   int64(a2)
       ),
       DBSNP
     ),
     DBSNP
   )")
  remove_versions('DBSNP')
}

ARRAYS[["VARIANT"]] = "
<rsid:               int64,
 ref:                string,
 alt:                string,
 site_quality:       float,
 filter:             string,
 info:               string, 
 original_signature: string>
[chrom    = 1:25,1,0,
 pos      = 0:*,10000000,0, 
 synthetic=0:999,1000,0]"

get_minimal_representation = function(pos, ref, alt)
{
  pos = as.numeric(pos)
  if(nchar(ref)==1 && nchar(alt)==1)
  {
    return(c(pos,ref,alt))
  }
  while( substr(alt, nchar(alt), nchar(alt)) == substr(ref, nchar(ref), nchar(ref)) && min(nchar(ref), nchar(alt)) > 1)
  {
    alt = substr(alt, 1, nchar(alt)-1)
    ref = substr(ref, 1, nchar(ref)-1)
  }
  while( substr(alt,1,1) == substr(ref,1,1) && min(nchar(ref), nchar(alt)) > 1)
  {
    alt = substr(alt, 2, nchar(alt))
    ref = substr(ref, 2, nchar(ref))
    pos  = pos +1
  }
  return(c(pos, ref, alt))
}


#Streaming parser for a block of variant data. This is run in DB in parallel. Primary objective is to pull apart ref/alt sequences.
#TODO: add pertinent fields from info if needed
#INPUT: chrom:string, pos:string, rsid: string, ref:string, alt:string, qual:string, filter:string, info:string
#OUTPUT: chrom:string, pos:string, rsid: string, ref:string, alt:string, qual:string, filter:string, original_signature:string, info:string
#Output is disassembled into minimal representations
#Test: 
# v= data.frame(chrom= c('1','2'), pos=c('123', '456'), rsid=c('rs1','rs2'), ref=c('AAT', 'A'), alt=c('A,T', 'T'), qual=c('100', '.'), filter=c('PASS', 'PASS'), info=c('INFOA','INFOB'), stringsAsFactors=FALSE)
# parse_variants(v)
parse_variants = function(variant_frame)
{
  in_chrom = variant_frame[, 1]
  in_pos = variant_frame[, 2]
  in_rsid = variant_frame[, 3]
  in_ref = variant_frame[, 4]
  in_alt = variant_frame[, 5]
  in_qual= variant_frame[, 6]
  in_filter=variant_frame[, 7]
  in_info = variant_frame[, 8]
  out_chrom = c()
  out_pos = c()
  out_rsid = c()
  out_ref = c()
  out_alt = c()
  out_qual = c()
  out_filter= c()
  out_info = c()
  out_original_signature = c()
  for(i in 1:nrow(variant_frame))
  {
    alternates = strsplit(in_alt[i], ',', fixed=T)[[1]]
    for(alt in alternates)
    {
      mr = get_minimal_representation(in_pos[i], in_ref[i], alt)
      out_chrom = c(out_chrom, in_chrom[i])
      out_pos   = c(out_pos, mr[1])
      out_rsid  = c(out_rsid, in_rsid[i])
      out_ref   = c(out_ref, mr[2])
      out_alt   = c(out_alt, mr[3])
      out_qual  = c(out_qual, in_qual[i])
      out_filter = c(out_filter, in_filter[i])
      out_info = c(out_info, in_info[i])
      sig = paste(in_chrom[i], in_pos[i], in_ref[i], in_alt[i], sep="-")
      out_original_signature = c(out_original_signature, sig)
    }
  }
  return(data.frame(out_chrom,out_pos,out_rsid,out_ref,out_alt,out_qual,out_filter,out_info, out_original_signature, stringsAsFactors = FALSE))
}

load_variants = function()
{
  connect()
  #TODO: make minimal representation
  #TODO: do we want to parse the info field here?
  run_pipe(SITES_VCF_FILE, path='/tmp/fifo', command='zcat')
  iquery(DB, "
   insert(
     redimension(
      apply(
        stream(
         filter(
          aio_input('/tmp/fifo', 'num_attributes=8'),
          substr(a0, 0,1)<>'#'
         ),
         'Rscript /home/scidb/GlobalBiobankEngine/scidb/streaming_parser.R parse_variants',
         'format=df',
         'types=string,string,string,string,string,string,string,string,string'
        ),
        chrom,  int64(a0),
        pos,    int64(a1),
        rsid,   iif(strlen(a2)>1, int64(substr(a2, 2, 100)), int64(null)),
        ref,    a3,
        alt,    a4, 
        site_quality,   dcast(a5, float(null)),
        filter, a6,
        info,   a7,
        original_signature, a8
       ),
       VARIANT
     ),
     VARIANT
    )")
  remove_versions('VARIANT')
}

ARRAYS[["ICD_INDEX"]] = "
<icd:    string>
[icd_id =0:*,1000000,0]"

map_icd_to_id = function(icd_ids)
{
  connect()
  filter_cond = paste0("icd='", icd_ids, "'", collapse=" or ")
  rs = iquery(DB, sprintf("filter(ICD_INDEX, %s)", filter_cond), return=T)
  existing_icds = rs$icd
  existing_ids = rs$icd_id
  if(nrow(rs) != length(icd_ids))
  {
    next_id = iquery(DB, "aggregate(apply(ICD_INDEX, i, icd_id), max(i) as icd_id)", return=T)$icd_id
    if (is.na(next_id))
    {
      next_id=0
    } else
    {
      next_id=next_id+1
    }
    for( icd in icd_ids)
    {
      if( !(icd %in% existing_icds))
      {
        print(paste("Adding ICD", icd, "ID", next_id))
        iquery(DB, sprintf("insert(redimension(build(<icd:string>[icd_id=%i:%i], '%s'), ICD_INDEX), ICD_INDEX)", next_id, next_id, icd))    
        existing_icds = c(existing_icds, icd)
        existing_ids = c(existing_ids, next_id)
        next_id = next_id+1
      }
    }
    remove_versions("ICD")
  }
  return(data.frame(existing_icds, existing_ids))
}

ARRAYS[["ICD"]] = "
<affyid:             string,
 or_val:             double,
 se:                 double,
 pvalue:             double,
 lor:                double,
 log10pvalue:        double,
 l95or:              double,
 u95or:              double>
[icd_id   = 0:*,20,0,
 chrom    = 1:25,1,0,
 pos      = 0:*,10000000,0,
 synthetic= 0:999,1000,0]"

load_icd = function()
{
  connect()
  n = 1 
  nfiles = length(ICD_STATS_FILES)
  filter1 = read.table(FILTER_FILE_1, header=T)
  filter1 = as.character(filter1[,1])
  filter2 = read.table(FILTER_FILE_2)
  filter2 = as.character(filter2[,1])
  filter = c(filter1, filter2)
  filter_scidb = as.scidb(DB, filter)
  while(n <= nfiles)
  {
    start = n
    end = min(n+7, length(ICD_STATS_FILES))
    files = ICD_STATS_FILES[start:end]
    n_loaded_files = length(files)
    icds = strsplit(basename(files), '.', fixed=T)
    icd_names = c()
    pipes = c()
    load_instances = c()
    for(i in 1:n_loaded_files)
    {
      icd_names = c(icd_names, icds[[i]][2])
      pipes = c(pipes, paste0("/tmp/fifo_",i))
      load_instances = c(load_instances, i-1)
    }
    icd_map = map_icd_to_id(icd_names)
    icd_ids = c()
    for(i in 1:n_loaded_files)
    {
      icd_ids = c(icd_ids, subset(icd_map, existing_icds==icd_names[i])$existing_ids)
    }
    icd_id_assignment = 'null'
    for(i in 1:n_loaded_files)
    {
      instance_id = load_instances[i]
      icd_id_assignment = sprintf("iif(src_instance_id = %i, %i, %s)", as.integer(instance_id), as.integer(icd_ids[i]), icd_id_assignment)
      run_pipe(files[i], path=pipes[i], command='zcat')
    }
    paths_string = paste0("'paths=",paste(pipes,collapse=";"), "'")
    inst_string = paste0("'instances=",paste(load_instances,collapse=";"), "'")
    #TODO: filter out inconsequential positions
    #TODO: add the log calculation (can also do that post-query)
    #run_pipe(file, path='/tmp/fifo', command='zcat')
    #icd = strsplit(basename(file), '.', fixed=T)[[1]][2]
    #icd_id = map_icd_to_id(icd)
    print(paste("Loading", files, "id", icd_ids))
    iquery(DB, sprintf("
      insert(
        redimension(
          apply(
            filter(
             index_lookup(
              aio_input(%s, %s, 'num_attributes=12') as INPUT,
              %s as FILTER, 
              INPUT.a1,
              present_in_filter
             ),
             substr(a0, 0,1)<>'#' and present_in_filter is null and a6 = 'ADD' and a11<>'NA' and dcast(a11, double(null)) <> 0
            ),
            icd_id, %s,
            chrom,  int64(a0),
            pos,    int64(a1),
            affyid,  a2,
            or_val,  dcast(a8, double(null)),
            se,      dcast(a9, double(null)),
            pvalue,  dcast(a11, double(null)),
            lor,     log(dcast(a8, double(null))),
            log10pvalue,  log10(dcast(a11, double(null))) * -1,
            l95or,        exp( log(dcast(a8, double(null))) - 1.96 * dcast(a9, double(null))),
            u95or,        exp( log(dcast(a8, double(null))) + 1.96 * dcast(a9, double(null)))
          ),
          ICD
        ),
        ICD
      )", 
      paths_string,  inst_string, filter_scidb@name, icd_id_assignment)
    )
    remove_versions('ICD')
    n = n+8
  }
}

ARRAYS[["AFFYID_INDEX"]] = "
<affyid:string, chrom:int64, pos:int64>[affyid_int=0:*,10000000,0]"

ARRAYS[["RSID_INDEX"]]
"<chrom:int64, pos:int64>[rsid = 0:*,10000000,0]"

make_extra_indeces = function()
{
  iquery(DB, "
    store(
     redimension(
      apply(
       grouped_aggregate(
        apply(
         ICD,
         chrom, chrom, pos, pos
        ),
        max(chrom) as chrom, max(pos) as pos,
        affyid
       ),
       affyid_int, int64(substr(affyid, 5, 100))
      ),
      AFFYID_INDEX
     ),
     AFFYID_INDEX
    )")
  remove_versions('AFFYID_INDEX')
  
  #Put all the RSIDs in one pile: why not
  #TODO: whoa looks like DBSNP and VARIANT are different assemblies in my case
  #Maybe I need an updated variant file
  iquery(DB, "
    store(
     redimension(
--      TODO: re-enable this when I get an updated variant file
--      equi_join(
--       grouped_aggregate(apply(VARIANT, chrom, chrom, pos, pos), max(chrom) as chrom, max(pos) as pos, rsid), 
--       grouped_aggregate(apply(DBSNP, chrom, chrom, pos, pos), max(chrom) as chrom, max(pos) as pos, rsid), 
--       'left_names=chrom,pos,rsid', 'right_names=chrom,pos,rsid', 'left_outer=1', 'right_outer=1'
--      ),
      grouped_aggregate(apply(VARIANT, chrom, chrom, pos, pos), max(chrom) as chrom, max(pos) as pos, rsid), 
      RSID_INDEX
     ),
     RSID_INDEX
    )")
}


ARRAYS[["ICD_INFO"]] = "
<case:int64, icd_name:string>[icd_id=0:*,1000000,0]"

load_icd_info = function()
{
  connect()
  iquery(DB, sprintf(
    "insert(
      redimension(
       apply(
        index_lookup(
         aio_input(
          '%s', 'num_attributes=6'
         ),
         ICD_INDEX,
         a0,
         icd_id
        ),
        case, dcast(a1, int64(null)),
        icd_name, a2
       ),
       ICD_INFO, false
      ),
      ICD_INFO
    )",
    ICD_INFO_FILE
  ))
}