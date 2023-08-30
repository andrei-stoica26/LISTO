#' @export

AddSynonyms <- function(synonyms, geneList){
  lapply(geneList, function(x) synonyms[[x]] <- geneList)
  return (synonyms)
}

CheckSynonyms <- function(seuratObj, synonyms, gene){
  message(str_c("Gene: ",gene))
  print(intersect(rownames(seuratObj),synonyms[[gene]]))
}

CheckMarkersInRownames <- function(seuratObj, markerNames, synonyms)
  return(invisible(lapply(setdiff(markerNames, rownames(seuratObj)),
                          function(x) CheckSynonyms(seuratObj, synonyms,x))))

synonymLists <- list(
  c("CD93","ECSM3","C1QR1","MXRA4"),
  c("CD44","HUTCH-I","HCELL","CSPG8","MC56","MDU2","MDU3","MIC4","IN"),
  c("DCLK1","DCDC3A","KIAA0369","DCAMKL1","DCLK"),
  c("LGR5","GPR49","GPR67","HG38","FEX"),
  c("ABCG2","BCRP","ABCP","MXR","EST157481","CD338","BCRP1","EC","UAQTL1","ABC15","GOUT1",
    "MXR-1","BMDP","MXR1","MRX"),
  c("CXCR4","LESTR","HM89","D2S201E","NPYY3R","HSY3RR","NPY3R","CD184","NPYR","LAP-3","NPYRL",
    "FB22","LCR1","WHIMS1","CXC-R4","CXCR-4","WHIMS","LAP3", "WHIM"),
  c("NES", "FLJ21841"),
  c("ABCB1","PGY1","CD243","GP170","ABC20","P-170","MDR1","CLCS","P-GP"),
  c("POU5F1","OCT3","OTF3","OCT-4","OTF-3","OCT4","OTF4"),
  c("ITGB4","CD104","GP150","JEB5A","JEB5B"),
  c("CD123","IL-3RA","IL3R","IL3RAY","IL3RX","IL3RY"),
  c("FUT4","FUC-TIV","FCT3A","CD15","SSEA-1","FUTIV"),
  c("HNF1A","LFB1","HNF1","TCF1","HNF-1A","MODY3","TCF-1","IDDM20","HNF4A"),
  c("CENPA", "CENP-A"),
  c("CDCA3","TOME-1","GRCC8","TOME1","C8","MODY3"),
  c("DDIAS","C11orf82","FLJ38838","FLJ25416","NOXIN"),
  c("FAM72B", "P17"),
  c("HASPIN", "GSG2"),
  c("EID3", "EP300", "NSMCE4B", "NSE4B","FLJ25832", "NS4EB", "EID-3"),
  c("MIRN221", "MIR221"),
  c("F2RL2", "PAR3", "PAR-3"),
  c("NR0B1", "DAX1", "AHCH", "AHC", "DSS","DAX-1","NROB1","SRXY2","AHX","GTD","HHG"),
  c("MIRN21", "MIR21"),
  c("HIST1H1B", "H1-5","H1F5","H1B","H1"),
  c("HIST2H2AB","H2AC21","H2AB"),
  c("TM4SF19","OCTM4"),
  c("NID2","NID-2"),
  c("HIST1H3B","H3C2"),
  c("CPEB2","CPE-BP2","HCPEB-2","CPEB-2"),
  c("FLJ10213","EBLN2","EBLN-2","EBLN-1"),
  c("RSPO3","THSD2","FLJ14440","PWTSR","CRISTIN1","HPWTSR"),
  c("HIST1H2BF","H2BC7"),
  c("PROM1", "CORD12", "PROML1", "AC133", "CD133", "RP41", "MCDR2", "STGD4", "MSTP061"),
  c("ALDH1A1", "RALDH1", "ALDH1", "PUMB1", "ALDH-E1", "RALDH", "ALDHII", "ALDC", "ALDH11", "HEL-9", "HEL12"),
  c("KIT", "SCFR", "CD117", "PBT", "EC", "MASTC"),
  c("SOX2", "MCOPS3", "ANOP3"),
  c("REG4", "RELP", "GISP", "REG-IV", "REG-4"),
  c("AGR2", "HAG-2", "AG2", "PDIA17", "XAG-2", "AG-2", "HPC8", "HEL-5", "GOB-4"),
  c("H1-4", "HIST1H1E", "H1F4", "RMNS", "H1E"),
  c("S100P", "MIG9", "S100E")
)


synonyms <- hash()
for (elem in synonymLists) synonyms <- AddSynonyms(synonyms, elem)


