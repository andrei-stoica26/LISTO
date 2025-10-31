#' Compute the overlaps of different pathways for two enrichment results
#'
#' This function compute the overlaps of different pathways for
#' two \code{enrichResult} objects.
#'
#' @param joinDF A data frame joining two \code{enrichResult} objects generated
#' with \code{joinER}
#' @param mtMethod Multiple testing correction method. Choices are
#' Bonferroni ('bf'), Benjamini-Hochberg('bh'), and the default
#' Benjamini-Yekutieli ('by').
#' @param ... Additional arguments passed to the multiple testing correction
#' function. See \code{?hammers::mtCorrectDF} for information.
#'
#' @return A data frame of pathways with significant gene overlaps.
#'
#' @examples
#' m1 <- genesER(c('AURKA', 'PTTG2', 'MKI67', 'RRM2', 'BUB1', 'KIF20C'),
#' 'human')
#' m2 <- genesER(c('AURKA', 'TOP2A', 'CENPF', 'BIRC5', 'BUB1', 'KIF20C'),
#' 'human')
#' df <- joinER(m1@result, m2@result)
#' res <- erDFPathwayOverlap(df)
#'
#' @export
#'
erDFPathwayOverlap <- function(df, mtMethod = c('by', 'bf', 'bh'), ...){
    mtMethod <- match.arg(mtMethod, c('by', 'bf', 'bh'))
    pvals <- vapply(df$Description, function(x) {
        rowDF <- df[df$Description==x, ]
        genes1 <- str_split(rowDF[['geneID_1']], '/')[[1]]
        genes2 <- str_split(rowDF[['geneID_2']], '/')[[1]]
        n <- as.integer(str_split(rowDF[['BgRatio']], '/')[[1]][1])
        return(pvalOverlap(genes1, genes2, n))
    }, numeric(1))
    res <- data.frame(term = names(pvals),
                      pval = pvals)
    res <- mtCorrectDF(res, mtMethod, ...)
    return(res)
}
