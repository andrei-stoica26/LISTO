#' @export

ReadCountOutput <- function(dir, fileName = "spliced") {
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/", fileName, ".mtx"))
  m <- Matrix::t(m)
  colnames(m) <- readLines(file(paste0(dir, "/", fileName, ".barcodes.txt")))
  rownames(m) <- readLines(file(paste0(dir, "/", fileName, ".genes.txt")))
  return(m)
}

GetRawCounts <- function(sampleName, countsName = "counts_filtered_", fileName = "spliced"){
  m <- ReadCountOutput(str_c(countsName, sampleName), fileName)
  m <- m[Matrix::rowSums(m) > 0,]
  tr2g <- read_tsv("t2g.txt", col_names = c("transcript", "gene", "gene_symbol")) %>%
    select(-transcript) %>%
    distinct()
  rownames(m) <- tr2g$gene_symbol[match(rownames(m), tr2g$gene)]
  m <- m[!is.na(rownames(m)), ]
  return(m)
}

GetRawCountsUnspliced <- function(sampleName)
  GetRawCounts(sampleName, "counts_unfiltered_", "unspliced")

SaveCounts <- function(datasetName, samples){
  dataFolder <- str_c(datasetName, "Data/")
  dir.create(dataFolder)
  resMats <- lapply(samples, GetRawCounts)
  invisible(lapply(1:length(samples), function(i)
    saveRDS(as(resMats[[i]], "dgCMatrix"), str_c(dataFolder, datasetName, "_Counts_", samples[[i]], ".rds"))))
}





