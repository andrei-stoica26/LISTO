#' @importFrom AnnotationDbi select
#' @importFrom clusterProfiler enrichGO enrichKEGG enrichWP
#' @importFrom DOSE setReadable
#' @importFrom org.Dr.eg.db org.Dr.eg.db
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @importFrom stringr str_split
#'
NULL

#' Perform enrichment analysis on a set of Entrez IDs
#'
#' This function performs enrichment analysis on a set of Entrez IDs.
#'
#' @param entrezIDs A character vector of gene Entrez IDs.
#' @inheritParams entrezGenes
#' @param funString Name of enrichment function from \code{clusterProfiler}.
#' Must be a character selected from 'enrichGO', 'enrichKEGG' and 'enrichWP'.
#'
#' @return Enrichment result.
#'
#' @keywords internal
#'
getEnrichmentResult <- function(entrezIDs,
                                species = c('human', 'mouse', 'zebrafish'),
                                funString = c('enrichGO',
                                              'enrichKEGG',
                                              'enrichWP')){
    species <- match.arg(species, c('human', 'mouse', 'zebrafish'))
    funString <- match.arg(funString, c('enrichGO', 'enrichKEGG', 'enrichWP'))
    df <- data.frame(database = c('org.Hs.eg.db', 'org.Mm.eg.db',
                                  'org.Dr.eg.db'),
                     code = c('hsa', 'dre', 'mmu'))
    rownames(df) <- c('human', 'mouse', 'zebrafish')
    db <- df[species, 1]
    code <- df[species, 2]
    return(switch(funString,
                  "enrichGO" = enrichGO(entrezIDs, OrgDb=db, ont="BP",
                                        readable=TRUE, pvalueCutoff=0),
                  "enrichWP" = setReadable(enrichWP(entrezIDs, 'Homo sapiens'),
                                           db, 'ENTREZID'),
                  "enrichKEGG" = setReadable(enrichKEGG(entrezIDs, code),
                                             db, 'ENTREZID'),
    ))
}

#' Convert gene symbols to Entrez IDs.
#'
#' This function converts gene symbols to Entrez IDs.
#'
#' @param genes A character vector of gene symbols.
#' @param species Species. Must be one of 'human', 'mouse' and 'zebrafish'.
#'
#' @return The Entrez IDs of the genes.
#'
#' @keywords internal
#'
entrezGenes <- function(genes, species = c('human', 'mouse', 'zebrafish')){
    species <- match.arg(species, c('human', 'mouse', 'zebrafish'))
    return(switch(species,
                  'human' = AnnotationDbi::select(org.Hs.eg.db,
                                                  keys=genes,
                                                  columns=c("ENTREZID", "SYMBOL"),
                                                  keytype="SYMBOL")[[2]],
                  'zebrafish' = AnnotationDbi::select(org.Dr.eg.db,
                                                      keys=genes,
                                                      columns=c("ENTREZID", "SYMBOL"),
                                                      keytype="SYMBOL")[[2]],
                  'mouse' = AnnotationDbi::select(org.Mm.eg.db, keys=genes,
                                                  columns = c("ENTREZID", "SYMBOL"),
                                                  keytype="SYMBOL")[[2]]))
}

#' Perform enrichment analysis on a set of genes
#'
#' This function performs enrichment analysis on a set of genes.
#'
#' @inheritParams entrezGenes
#' @inheritParams getEnrichmentResult
#' @param pvalThr p-value threshold.
#' @param addNegLog Whether to compute the negative log10 of the adjusted
#' p-value.
#' @param pvalOffset Offset used to avoid zeros inside the logarithm function.
#' Ignored if \code{addNegLog} is \code{FALSE} (as default).
#'
#' @return Enrichment result.
#'
#' @examples
#' m <- genesER(c('AURKA', 'TOP2A', 'CENPF', 'PTTG2', 'MKI67', 'BIRC5', 'RRM2'),
#' 'human')
#'
#' @export
#'
genesER <- function(genes,
                    species,
                    funString = c('enrichGO','enrichKEGG', 'enrichWP'),
                    pvalThr = 0.05,
                    addNegLog = FALSE,
                    pvalOffset = 1e-317
){
    m <- getEnrichmentResult(entrezGenes(genes, species), species, funString)
    m@result <- m@result[m@result$p.adjust < pvalThr, ]
    if(addNegLog)
        m@result$nlog.padj <- -log(m@result$p.adjust + pvalOffset, base=10)
    return(m)
}

#' Extract genes enriched for terms
#'
#' This function extracts genes enriched for term from an \code{enrichResult}
#' object.
#'
#' @param er Enrichment result of the \code{enrichResult} class.
#' @param terms Terms for which enriched genes should be extracted.
#'
#' @return A character vector of genes enriched for the input terms.
#'
#' @keywords internal
#'
termGenesHelper <- function(er, terms)
    return(unique(unlist(lapply(er@result[er@result$Description %in%
                                              terms, ]$geneID,
                                function(x) str_split(x, "/")[[1]]))))

#' Extract genes enriched for terms
#'
#' This function extracts genes enriched for terms from an \code{enrichResult}
#' object.
#'
#' @inheritParams termGenesHelper
#' @param negTerms Terms for which enriched genes should be subtracted from the
#' genes enriched for \code{terms}.
#'
#' @return Genes enriched for terms.
#'
#' @examples
#' m <- genesER(c('AURKA', 'TOP2A', 'CENPF', 'PTTG2', 'MKI67', 'BIRC5', 'RRM2'),
#' 'human')
#' termGenes(m, 'chromosome segregation', 'meiosis I')
#'
#' @export
#'
termGenes <- function(er, terms, negTerms = NULL){
    posGenes <- termGenesHelper(er, terms)
    if(!is.null(negTerms))
        negGenes <- termGenesHelper(er, negTerms) else
            negGenes <- NULL
        return(sort(setdiff(posGenes, negGenes)))
}

#' Join the data frames of two enrichment results
#'
#' This function joins the data frames of two \code{enrichResult} objects.
#'
#' @param df1 The \code{result} data frame of an \code{enrichResult} object.
#' @param df2 The \code{result} data frame of an \code{enrichResult} object.
#' @param sharedCols Columns from the \code{enrichResult} objects used for the
#' join that have the same values for each term shared by the two objects.
#' @param specificCols Columns from the \code{enrichResult} objects used for the
#' join that have potentially different values for terms shared by the
#' two objects.
#'
#' @return Genes enriched for terms.
#'
#' @examples
#' m1 <- genesER(c('AURKA', 'PTTG2', 'MKI67', 'RRM2', 'BUB1', 'KIF20C'),
#' 'human')
#' m2 <- genesER(c('AURKA', 'TOP2A', 'CENPF', 'BIRC5', 'BUB1', 'KIF20C'),
#' 'human')
#' df <- joinER(m1@result, m2@result)
#'
#'
#' @export
#'
joinER <- function(df1,
                   df2,
                   sharedCols = c('Description', 'BgRatio'),
                   specificCols = c('p.adjust', 'geneID')){
    if (length(setdiff(sharedCols, c('ID', 'Description', 'BgRatio'))))
        stop('`sharedCols` must be included in ',
             'c("ID", "Description", "BgRatio")')
    if (length(setdiff(specificCols, c('GeneRatio', 'pvalue', 'p.adjust',
                                       'qvalue', 'geneID', 'Count'))))
        stop('`specificCols` must be included in ',
             'c("GeneRatio", "pvalue", "p.adjust", ',
             '"qvalue", "geneID", "Count")')

    sharedIDs <- intersect(rownames(df1), rownames(df2))

    dfTemp1 <- df1[sharedIDs, specificCols]
    dfTemp2 <- df2[sharedIDs, specificCols]
    colnames(dfTemp1) <- paste0(specificCols, '_1')
    colnames(dfTemp2) <- paste0(specificCols, '_2')
    res <- do.call(cbind, list(df1[sharedIDs, sharedCols],
                               dfTemp1,
                               dfTemp2))

    return(res)
}
