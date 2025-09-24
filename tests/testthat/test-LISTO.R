test_that("subsetOverlapProbMNk", {
    res <- sum(vapply(seq(0, 70), function(i)
        pOverlapMNk(300, 70, 110, i), numeric(1)))
    expect_equal(res, 1, tolerance=0.0001)
})
