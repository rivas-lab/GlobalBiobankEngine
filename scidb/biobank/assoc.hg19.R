source("/usr/local/src/biobank/pkg/R/biobank.R")

bb <- connect(username = "scidbadmin", password = "Paradigm4")
ns_name <- "RIVAS_HG19"
as_name <- "UK_Biobank_Array_European"
skip_list <- c()
#  "ukb24983_v2_hg19.HC383.genotyped.PHENO1.glm.logistic.hybrid.gz")
skip_filename <- "skip.list"

# If a skip file does not exist,
# assume we need to create the association set

if (file.exists(skip_filename)) {
  skip_list <- c(skip_list, readLines(skip_filename))
} else {
#   bb$delete_association_set(ns_name, as_name)
#  bb$create_association_set(ns_name, as_name, format(Sys.time(), "%Y-%m-%d"))
}

directory <- "/opt/biobankengine/GlobalBioBankEngineRepo/gbe_data/20200717_plink/"
pattern <- paste0(
  "([a-zA-Z0-9_]+)",
  ".metal.plink.tsv.gz")

for (filename in list.files(path = directory, pattern = pattern))
{
  phenotype <- gsub(
    pattern,
    "\\1",
    filename)

  print(paste("Processing", filename, "phenotype", phenotype))

  if (filename %in% skip_list) {
    print(paste0("Skipping ",filename))
    next
  }

  if (grepl("INI", filename) | grepl("QT", filename)) {
    col.names = c(
      "chrom",
      "pos",
      "rsid",        # ID
      "ref",
      "alt",
      "a1",
      "test",
      "nobs",        # OBS_CT
      "beta",
      "se",
      "t_stat",
      "pvalue")
    colClasses = c(
      "integer",
      "integer",
      "character",
      "character",
      "character",
      "character",
      "character",
      "integer",
      "numeric",
      "numeric",
      "numeric",
      "numeric")
  } else { # ".logistic."
    col.names = c(
      "chrom",
      "pos",
      "rsid",        # ID
      "ref",
      "alt",
      "a1",
      "firth",
      "test",
      "nobs",        # OBS_CT
      "odds_ratio",  # OR
      "se",
      "z_stat",
      "pvalue")
    colClasses = c(
      "integer",
      "integer",
      "character",
      "character",
      "character",
      "character",
      "character",
      "character",
      "integer",
      "numeric",
      "numeric",
      "numeric",
      "numeric")
  }

  #data <- read.delim(
  #  gzfile(paste0(directory, "/", filename)),
  #  header = FALSE,
  #  comment.char = "#",
  #  col.names = col.names,
  #  colClasses = colClasses)

  t2=proc.time()
  data  <- data.table::fread(
     cmd = paste0("zcat ", directory, "/", filename),
     header = TRUE,
     col.names = col.names,
     colClasses = colClasses);
  setattr(data, 'class', 'data.frame')
  if(class(data$pvalue) != 'numeric')
  {
     data$pvalue = as.numeric(data$pvalue)
  }
  print("Read finished")
  print(proc.time()-t2)


  if (grepl("INI", filename) | grepl("QT", filename)) {
    data$odds_ratio <- as.numeric(NA)
    data$z_stat <- as.numeric(NA)
    data$note <- ""
  } else { # ".logistic."
    data$beta <- as.numeric(NA)
    data$t_stat <- as.numeric(NA)
    data$note <- paste0("FIRTH=", data$firth)
    data$firth <- NULL
  }
  data$phenotype <- phenotype

  data$rsid <- NULL
  data$a1 <- NULL
  data$test <- NULL
data = subset(data, pvalue < 0.1)
  bb$upload_association_data(
    ns_name,
    as_name,
    data,
    allow_new_phenotypes = TRUE,
    allow_new_variants = TRUE,
    verbose = TRUE)
  gc()

  cat(paste("Done processing", filename, "\n"))
}

print("Computing thresholds")
bb$compute_pvalue_thresholds(ns_name, as_name)

# TODO Move out of the association set loader
print("Create RSID index")
bb$populate_rsid_index(ns_name)

print("Complete")
