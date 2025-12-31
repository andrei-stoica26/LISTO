# A self-contained script to generate the data in extdata
#
createExtData <- function(){
    if (requireNamespace(c('LISTO', 'qs2', 'scRNAseq','scuttle', 'Seurat',
                           'withr'), quietly=TRUE)){
        scObj <- scRNAseq::BaronPancreasData('human')
        scObj <- scuttle::logNormCounts(scObj)
        scObj <- Seurat::as.Seurat(scObj)

        donorMarkers <- LISTO::buildSeuratMarkerList(scObj, 'donor')
        labelMarkers <- LISTO::buildSeuratMarkerList(scObj, 'label')
        acinarMarkers <- c('PRSS1', 'KLK1', 'CTRC', 'PNLIP', 'AKR1C3', 'CTRB1',
                           'DUOXA2', 'ALDOB', 'REG3A', 'SERPINA3', 'PRSS3',
                           'REG1B','CFB',	'GDF15', 'MUC1','ANPEP', 'ANGPTL4',
                           'OLFM4','GSTA1', 'LGALS2', 'PDZK1IP1', 'RARRES2',
                           'CXCL17','UBD', 'GSTA2', 'LYZ', 'RBPJL', 'PTF1A',
                           'CELA3A','SPINK1', 'ZG16', 'CEL', 'CELA2A', 'CPB1',
                           'CELA1','PNLIPRP1', 'RNASE1', 'AMY2B', 'CPA2','CPA1',
                           'CELA3B','CTRB2', 'PLA2G1B', 'PRSS2', 'CLPS',
                           'REG1A', 'SYCN')

        alphaMarkers <- c('GCG', 'TTR', 'PCSK2', 'FXYD5', 'LDB2', 'MAFB',
                          'CHGA', 'SCGB2A1', 'GLS', 'FAP', 'DPP4', 'GPR119',
                          'PAX6', 'NEUROD1', 'LOXL4', 'PLCE1', 'GC', 'KLHL41',
                          'FEV', 'PTGER3', 'RFX6', 'SMARCA1', 'PGR', 'IRX1',
                          'UCP2', 'RGS4', 'KCNK16', 'GLP1R', 'ARX', 'POU3F4',
                          'RESP18', 'PYY', 'SLC38A5', 'TM4SF4', 'CRYBA2',
                          'SH3GL2','PCSK1', 'PRRG2', 'IRX2', 'ALDH1A1','PEMT',
                          'SMIM24','F10', 'SCGN', 'SLC30A8')

        betaMarkers <- c('INS', 'IAPP', 'GJD2', 'PDX1', 'SLC2A2', 'NPY', 'MAFA',
                         'PFKFB2','HOPX', 'PAX6', 'MAFB', 'CASR', 'EDARADD',
                         'SCGB2A1','TGFBR3','ADCYAP1', 'SH3GL2', 'NEUROD1',
                         'ISL1', 'RGS16','SMAD9', 'SIX3','BMP5', 'PIR',
                         'STXBP5', 'DLK1', 'MEG3','GCGR', 'LMX1A', 'JPH3',
                         'CD40', 'HAMP','EZH1', 'NTRK1', 'FXYD2', 'RIMS1',
                         'EFNA5', 'NPTX2', 'PAX4', 'PCSK2', 'G6PC2', 'SLC30A8',
                         'PCSK1', 'SCGN', 'IGF2', 'SYT13', 'FFAR2', 'SIX2')

        deltaMarkers <- c('SST', 'GABRB3', 'FRZB', 'MS4A8', 'BAIAP3', 'CASR',
                          'BCHE', 'UNC5B','EDN3', 'GHSR', 'PCSK1', 'GABRG2',
                          'POU3F1', 'BHLHE41', 'EHF', 'LCORL', 'ETV1', 'PDX1',
                          'LEPR', 'UCP2', 'NPTX2', 'FXYD2', 'IAPP', 'KCNK16',
                          'SCGN', 'ISL1', 'HHEX', 'RESP18', 'PAX4', 'RBP4',
                          'PCSK9', 'FFAR4')

        ductalMarkers <-  c('CFTR', 'SERPINA5', 'SLPI', 'TFF1', 'CFB',
                            'LGALS4', 'CTSH',	'PERP', 'PDLIM3', 'WFDC2',
                            'SLC3A1', 'AQP1', 'ALDH1A3', 'VTCN1',	'KRT19',
                            'TFF2', 'KRT7', 'CLDN4', 'LAMB3', 'TACSTD2',
                            'CCL2', 'DCDC2','CXCL2', 'CLDN10', 'HNF1B',
                            'KRT20', 'MUC1', 'ONECUT1', 'AMBP', 'HHEX',
                            'ANXA4', 'SPP1', 'PDX1', 'SERPINA3', 'GDF15',
                            'AKR1C3', 'MMP7', 'DEFB1', 'SERPING1', 'TSPAN8',
                            'CLDN1', 'S100A10', 'PIGR')

        gammaMarkers <- c('PPY', 'ABCC9', 'FGB', 'ZNF503', 'MEIS1', 'LMO3',
                          'EGR3', 'CHN2', 'PTGFR', 'ENTPD2', 'AQP3', 'THSD7A',
                          'CARTPT', 'ISL1', 'PAX6', 'NEUROD1', 'APOBEC2',
                          'SEMA3E', 'SLITRK6','SERTM1', 'PXK', 'PPY2P', 'ETV1',
                          'ARX', 'CMTM8', 'SCGB2A1', 'FXYD2', 'SCGN')
        signatures <- list(acinarMarkers, alphaMarkers, betaMarkers,
                           deltaMarkers, ductalMarkers, gammaMarkers)

        qs2::qs_save(donorMarkers, 'inst/extdata/donorMarkers.qs2')
        qs2::qs_save(labelMarkers, 'inst/extdata/labelMarkers.qs2')
        qs2::qs_save(signatures, 'inst/extdata/signatures.qs2')
        qs2::qs_save(rownames(scObj), 'inst/extdata/universe1.qs2')
        qs2::qs_save(Reduce(union, signatures), 'inst/extdata/universe2.qs2')

        seuratObj <- withr::with_seed(1, scuttle::mockSCE(ngenes=500))
        seuratObj <- Seurat::as.Seurat(seuratObj, data=NULL)
        seuratObj <- Seurat::NormalizeData(seuratObj)
        qs2::qs_save(seuratObj, 'inst/extdata/seuratObj.qs2')
    }
}

createExtData()
