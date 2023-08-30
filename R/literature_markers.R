#' @export

AddGeneSets <- function(seuratObj, tradMarkers, LitMarkers, sideLogMarkers)
  return(list(
    trad = intersect(tradMarkers, rownames(seuratObj)),
    Lit = lapply(LitMarkers, function(x) intersect(x, rownames(seuratObj))),
    side = sideLogMarkers[intersect(rownames(sideLogMarkers), rownames(seuratObj)),,drop = F]
  ))
tradMarkersPatient <- c("MET", "CD24", "CD44", "CD93", "PROM1", "DCLK1", "EPCAM","NANOG", "LGR5", "ABCG2", "ALDH1A1",
                        "LAP3",  "NES", "ABCB1", "POU5F1", "THY1", "KIT","SOX2", "TSPAN8","REG4",
                        "CD9", "EZH2", "FUT4", "SSEA-4", "AGR2", "GLRX3", "HNF4A")

tradMarkersA13A <- c("MET", "CD24", "CD44", "CD93", "PROM1", "DCLK1", "EPCAM", "NANOG", "LGR5", "ABCG2", "ALDH1A1",
                     "LAP3",  "NES", "ABCB1", "POU5F1", "THY1", "KIT", "SOX2", "TSPAN8","REG4",
                     "CD9", "EZH2", "FUT4", "SSEA-4", "AGR2", "GLRX3", "HNF4A")

prognosis <- c("ANLN", "ZWINT", "CEP55","TOP2A", "UBE2C", "ZWILCH", "CDK1", "STIL", "PCLAF",
               "GINS1", "CENPF", "PRC1", "RRM2", "ASPM", "FANCI", "KIF20A",
               "CENPU", "NUSAP1", "CENPK", "TK1")

gastric <- c("BUB1", "BUB1B", "NCAPH", "KIF14", "RACGAP1", "RAD54L", "TPX2", "KIF15", "KIF18B", "CENPF", "TTK",
             "KIF4A", "SGO2", "PLK4", "XRCC2", "C1orf112")

lung <- c("BUB1B", "CDC25A", "CDCA5", "CENPA", "DKC1", "NCAPH", "SKA3", "SPAG5","TIMELESS", "RAD51")

ovarian <- c("BUB1", "CDC20", "CCNB2", "DLGAP5", "KIF4A", "NEK2", "NUSAP1")

breast <- c("TPX2", "HJURP", "CDCA8", "PLK1", "KIFC1", "CENPA", "CCNB2", "KIF2C", "EXO1", "TTK", "KIF4A", "CDC25A",
            "MELK", "NDC80", "NCAPG", "CEP55", "NCAPH", "RAD54L", "KIF20A", "KIF18B", "ORC1", "CDC45", "KIF23",
            "CDC20","BUB1", "AURKB", "SKA1", "FOXM1", "SGO1", "DLGAP5", "CDCA3", "BUB1B")

liver <- c("KIF4A", "TTK", "CCNB1", "CDC20", "NCAPG", "CCNB2", "CDC45", "UBE2C", "CENPA", "AURKB", "RRM2", "CDCA8",
           "BIRC5", "TPX2","KIF2C")

glioma <- c("GINS2", "DBF4", "OIP5", "PBK", "CCNB2", "C1orf112", "CDCA8", "MELK", "DEPDC1B", "DDIAS", "RCC1",
            "DLGAP5", "NUF2", "AC099850.3", "FAM72B", "RAD51", "ORC1", "CDK1", "SGO1", "CDCA2", "KIF2C",
            "SKA3","PRR11", "BUB1", "TTK", "ESCO2", "FBXO5", "NCAPG","NDC80", "SGO2", "KIF4A", "TPX2",
            "NCAPH", "CKAP2L", "HASPIN", "CENPI", "SMC4", "CKAP2", "KIF15", "TOP2A", "ZWINT", "KIFC1", "LMNB1",
            "E2F2", "KNL1", "KIF14", "ASPM", "C18orf54", "MCM10", "KIF11", "CENPF")

pancreas <- c("NEK2", "PBK", "NCAPH", "CENPA", "TPX2", "PLK1", "CDC20",
              "KIF4A", "MKI67", "HJURP", "CKS2", "CCNA2", "KIF11", "ZWINT", "DTL", "UBE2C", "CDCA5",
              "GINS1", "CDKN3", "PTTG1", "RAD51", "CCNB2", "CDK1", "GINS2", "KIFC1", "SKA3", "NUF2",
              "CEP55", "BUB1", "KIF18B", "CDC45", "BIRC5", "ASF1B", "AURKA", "E2F1", "UBE2T")

signature <- c("TTK", "CDC20", "PBK", "NDC80", "KIF4A", "MELK", "PRC1", "KIF20A", "ECT2", "DTL", "KIF2C",
               "GPSM2", "OIP5", "PCLAF")

endometrial <- c("ORC6", "C1orf112", "RAD54L", "SGO2", "BUB1", "PLK4", "KIF18B", "BUB1B", "TTK", "NCAPG", "XRCC2",
                 "CENPF", "KIF15", "RACGAP1", "ARHGAP11A", "TPX2", "KIF14", "KIF4A", "NCAPH", "PLK1", "CDK1",
                 "MAD2L1")

bladder <- c("AURKA", "BUB1B", "CDCA5", "CDCA8", "KIF11", "KIF18B", "KIF2C", "KIFC1", "KPNA2", "NCAPG", "NEK2",
             "NUSAP1", "RACGAP1")

colon <- c("CHEK1", "BUB1", "KIF18A", "TTK", "PLK4", "NUP107", "SPC25", "DNA2", "DDIAS", "MCM10", "RFC4", "NCAPG",
           "BUB1B", "SUV39H2", "NCAPH", "KIF23", "CDK1", "MELK", "DEPDC1B", "NEIL3", "MTFR2", "PNPT1", "ORC6",
           "CCNA2", "MAD2L1", "CENPA", "XRCC2")

biomarkers <- c("KIF4A", "RRM2", "FAM83D", "ASPM", "NCAPG", "TPX2", "NUSAP1", "RACGAP1", "MKI67", "KIF20A", "CENPF",
                "UBE2C", "BUB1", "KIF11")

LitMarkers <- list(signature, pancreas, gastric, lung, ovarian, breast, glioma, endometrial, bladder,
                   colon, liver, prognosis, biomarkers)

allLit <- Reduce(union, LitMarkers)
LitNames <- list("Signature", "Pancreas", "Gastric", "Lung", "Ovarian", "Breast", "Glioma", "Endometrial", "Bladder",
                 "Colon", "Liver", "Prognosis", "Biomarkers", "Union")
LitMarkers <- append(LitMarkers, list(allLit))

sideLogMarkers <- data.frame(avg_log2FC = c(3.787, 2.48, 1.771, 1.742, 1.729, 1.589, 1.439, 1.436, 1.427, 1.377,
                                            1.349, 1.338, 1.321, 1.321, 1.307,
                                            1.298, 1.295, 1.284, 1.276, 1.217, 1.216, 1.152, 1.14, 1.137, 1.123,
                                            1.11, 1.104, 1.104, 1.095, 1.087,
                                            1.083, 1.082, 1.077, 1.074, 1.06, 1.038, 1.027, 1.025, 1.013))



rownames(sideLogMarkers) <- c("AKR1B10", "ABCG2", "EP300", "MIRN221", "F2RL2", "GDF15", "NR0B1", "GPAT3", "MIRN21", "MAP1B",
                              "CYP4F3", "HIST1H1B", "LRRFIP1", "SNORD122", "KIAA0319", "H2AC21", "H1-4", "TXNRD1", "TM4SF19",
                              "NID2", "NTS", "WNT5A", "SNORA1", "MAN1A1", "S100P","SNORA16B", "AKR1C3", "GDA", "HIST1H3B",
                              "LIMCH1", "SNORA25", "SNORD53", "CPEB2", "LIPH", "FLJ10213", "SNORD5", "RSPO3", "TMEM156",
                              "HIST1H2BF")

#For patient:

rownames(sideLogMarkers) <- c("AKR1B10", "ABCG2", "EP300", "MIRN221", "F2RL2", "GDF15", "NR0B1", "GPAT3", "MIRN21", "MAP1B",
                              "CYP4F3", "HIST1H1B", "LRRFIP1", "SNORD122", "KIAA0319", "HIST2H2AB", "HIST1H1E", "TXNRD1", "TM4SF19",
                              "NID2", "NTS", "WNT5A", "SNORA1", "MAN1A1", "S100P","SNORA16B", "AKR1C3", "GDA", "HIST1H3B",
                              "LIMCH1", "SNORA25", "SNORD53", "CPEB2", "LIPH", "EBLN2", "SNORD5", "RSPO3", "TMEM156",
                              "HIST1H2BF")

sideLogMarkersA13A <- data.frame(avg_log2FC = c(3.787, 2.48, 1.771, 1.742, 1.729, 1.589, 1.439, 1.436, 1.427, 1.377,
                                                1.349, 1.338, 1.321, 1.321, 1.307,
                                                1.298, 1.295, 1.284, 1.276, 1.217, 1.216, 1.152, 1.14, 1.137, 1.123,
                                                1.11, 1.104, 1.104, 1.095, 1.087,
                                                1.083, 1.082, 1.077, 1.074, 1.06, 1.038, 1.027, 1.025, 1.013))

rownames(sideLogMarkersA13A) <- c("AKR1B10", "ABCG2", "EP300", "MIRN221", "F2RL2", "GDF15", "NR0B1", "GPAT3", "MIRN21", "MAP1B",
                                  "CYP4F3", "HIST1H1B", "LRRFIP1", "SNORD122", "KIAA0319", "H2AC21", "H1-4", "TXNRD1", "TM4SF19",
                                  "NID2", "NTS", "WNT5A", "SNORA1", "MAN1A1", "S100P","SNORA16B", "AKR1C3", "GDA", "HIST1H3B",
                                  "LIMCH1", "SNORA25", "SNORD53", "CPEB2", "LIPH", "FLJ10213", "SNORD5", "RSPO3", "TMEM156",
                                  "HIST1H2BF")

sideMarkersA13A <- c("AKR1B10", "ABCG2", "EP300", "MIRN221", "F2RL2", "GDF15", "NR0B1", "GPAT3", "MIRN21", "MAP1B",
                     "CYP4F3", "HIST1H1B", "LRRFIP1", "SNORD122", "KIAA0319", "H2AC21", "H1-4", "TXNRD1", "TM4SF19",
                     "NID2", "NTS", "WNT5A", "SNORA1", "MAN1A1", "S100P","SNORA16B", "AKR1C3", "GDA", "HIST1H3B",
                     "LIMCH1", "SNORA25", "SNORD53", "CPEB2", "LIPH", "FLJ10213", "SNORD5", "RSPO3", "TMEM156",
                     "HIST1H2BF")
sideMarkersPatient <- c("AKR1B10", "ABCG2", "EP300", "MIRN221", "F2RL2", "GDF15", "NR0B1", "GPAT3", "MIRN21", "MAP1B",
                        "CYP4F3", "HIST1H1B", "LRRFIP1", "SNORD122", "KIAA0319", "HIST2H2AB", "HIST1H1E", "TXNRD1", "TM4SF19",
                        "NID2", "NTS", "WNT5A", "SNORA1", "MAN1A1", "S100P","SNORA16B", "AKR1C3", "GDA", "HIST1H3B",
                        "LIMCH1", "SNORA25", "SNORD53", "CPEB2", "LIPH", "EBLN2", "SNORD5", "RSPO3", "TMEM156",
                        "HIST1H2BF")


