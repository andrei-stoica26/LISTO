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

test_that("buildSeuratMarkerList works", {
    res <- buildSeuratMarkerList(seuratObj, 'Cell_Cycle', logFCThr=0.1)
    expect_true(is(res$G0, 'data.frame'))
    expect_equal(res$G0$p_val_adj[1], 0.8366238, tolerance=0.0001)
    expect_equal(length(res), 4)
})

