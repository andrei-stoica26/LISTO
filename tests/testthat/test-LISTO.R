test_that("checkNumColAll works", {
    df1 <- data.frame(fruit = c('apple', 'banana', 'cherry'),
                      cost = c(6, 5, 3))
    df2 <- data.frame(fruit = c('apple', 'banana', 'cherry'),
                      cost = c('6', '5', '3'))

    expect_no_message(checkNumCol(df1, 'cost'))
    expect_error(checkNumCol(df1, 'price'))
    expect_error(checkNumCol(df2, 'cost'))
    expect_no_message(checkNumCol(LETTERS, 'cost'))

    expect_no_message(checkNumColAll (list(df1, LETTERS), 'cost'))
    expect_error(checkNumColAll (list(df1, LETTERS), 'price'))
    expect_error(checkNumColAll (list(df2, LETTERS), 'cost'))
    expect_no_message(checkNumColAll (list(LETTERS, LETTERS), 'cost'))
})

test_that("getObjectValues works", {
    df <- data.frame(fruit = c('apple', 'banana', 'cherry'),
                      cost = c(6, 5, 3))
    expect_identical(getObjectValues(df, 'cost'), c(6, 5, 3))
    expect_identical(getObjectValues(df, NULL), -Inf)
    expect_identical(getObjectValues(LETTERS, 'cost', FALSE), Inf)
})

test_that("generateCutoffs works", {
    res <- generateCutoffs(donorMarkers[[1]],
                           donorMarkers[[2]],
                           numCol='avg_log2FC')
    expect_equal(res[1], 2.743160, tolerance=0.0001)
    res <- generateCutoffs(donorMarkers[[2]],
                           donorMarkers[[3]],
                           labelMarkers[[1]],
                           numCol='avg_log2FC')
    expect_equal(res[1], 3.501303, tolerance=0.0001)

})

test_that("factorization functions work", {
    expect_identical(factorialPrimePowers(8), c(7, 2, 1, 1))
    expect_identical(powerProduct(c(2, 3, 5), c(4, 2, 6)), 2250000)
})

test_that("multiple testing functions work", {
    pvals <- c(0.032, 0.001, 0.0045, 0.051, 0.048)
    res <- mtCorrectV(pvals)
    expect_equal(res, c(0.11645000, 0.01141667, 0.02568750, 0.11645000,
                        0.11645000), tolerance=0.0001)
    res <- mtCorrectV(pvals, 'hochberg', 'median')
    expect_equal(res, 0.051, tolerance=0.0001)
    df <- data.frame(elem = c('A', 'B', 'C', 'D', 'E'),
    pval = pvals)
    res <- mtCorrectDF(df)
    expect_equal(res$pvalAdj, c(0.01141667, 0.02568750), tolerance=0.0001)
})

test_that("probCounts2MN works", {
    res <- probCounts2MN(300, 70, 110, 6)
    expect_equal(res, 2.091706e-09, tolerance=0.0001)
    res <- sum(vapply(seq(0, 70), function(i)
        probCounts2MN(300, 70, 110, i), numeric(1)))
    expect_equal(res, 1, tolerance=0.0001)
})

test_that("probCounts3N works", {
    res <- probCounts3N(25, 20, 30, 140, 8)
    expect_equal(res, 6.875155e-08, tolerance=0.0001)
    res <- sum(vapply(seq(0, 20), function(k) probCounts3N(25, 20, 30, 140, k),
                      numeric(1)))
    expect_equal(res, 1, tolerance=0.0001)
})

test_that("pvalCounts functions work", {
    expect_equal(pvalCounts2MN(300, 23, 24, 6), 0.005571074,
                 tolerance=0.0001)
    expect_equal(pvalCounts3N(300, 200, 250, 400, 180), 5.101079e-62,
                               tolerance=0.0001)
})

test_that("pvalObjects works", {
    res <- pvalObjects(LETTERS[seq(2, 7)],
                       LETTERS[seq(3, 19)],
                       universe1=LETTERS)
    expect_equal(res, 0.2956522, tolerance=0.0001)
    res <- pvalObjects(LETTERS[seq(2, 7)],
                       LETTERS[seq(3, 19)],
                       LETTERS[seq(4, 8)],
                       universe1=LETTERS,
                       type='3N')
    expect_equal(res, 0.0007643267, tolerance=0.0001)
    res <- pvalObjects(LETTERS[seq(2, 7)],
                       LETTERS[seq(3, 8)],
                       universe1=LETTERS[seq(1, 16)],
                       universe2=LETTERS[seq(2, 26)],
                       type='2MN')
    expect_equal(res, 0.01098901, tolerance=0.0001)
})

test_that("pvalSets functions work", {
    res <- pvalSets2N(LETTERS[seq(4, 10)], LETTERS[seq(7, 15)], LETTERS)
    expect_equal(res, 0.1585284, tolerance=0.0001)
    res <- pvalSets2MN(LETTERS[seq(4, 10)],
                       LETTERS[seq(7, 15)],
                       LETTERS[seq(19)],
                       LETTERS[seq(6, 26)])
    expect_equal(res, 0.3776224, tolerance=0.0001)
    res <- pvalSets3N(LETTERS[seq(4, 10)],
                      LETTERS[seq(7, 15)],
                      LETTERS[seq(19)],
                      LETTERS)
    expect_equal(res, 0.05096209, tolerance=0.0001)
})

test_that("runLISTO works", {
    res <-  runLISTO(donorMarkers, labelMarkers, universe1=universe1,
                     numCol='avg_log2FC', verbose=FALSE)
    expect_equal(dim(res), c(length(donorMarkers) * length(labelMarkers), 4))
    expect_equal(res[1, 4], 1.636216e-51, tolerance=0.0001)
    res <- runLISTO(donorMarkers, signatures, universe1=universe1,
                    universe2=universe2, numCol='avg_log2FC',
                    verbose=FALSE)
    expect_equal(dim(res), c(length(donorMarkers) * length(signatures), 4))
    expect_equal(res[1, 4], 9.317070e-41, tolerance=0.0001)
    res <-  runLISTO(donorMarkers, labelMarkers, signatures,
                     universe1=universe1, numCol='avg_log2FC',
                     verbose=FALSE)
    expect_equal(dim(res), c(length(donorMarkers) *
                                 length(labelMarkers) *
                                 length(signatures), 5))
    expect_equal(res[1, 4], 1.259111e-86, tolerance=0.0001)
})

test_that("vector functions work", {
    res <- vSum(c(1, 4), c(2, 8, 6), c(1, 7), c(10, 4, 6, 7))
    expect_identical(res, c(14, 23, 12, 7))
    res <- vChoose(8, 4)
    expect_identical(res, c(1, 0, 1, 1))
    res <- vNumeratorMN(20, 8, 6, 2)
    expect_identical(res, c(2, 2, 1, 1, 1))
})

test_that("buildSeuratMarkerList works", {
    res <- buildSeuratMarkerList(seuratObj, 'Cell_Cycle', logFCThr=0.1)
    expect_true(is(res$G0, 'data.frame'))
    expect_equal(res$G0$p_val_adj[1], 0.8366238, tolerance=0.0001)
    expect_equal(length(res), 4)
})



