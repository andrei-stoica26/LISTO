test_that("probCounts2MN works", {
    res <- sum(vapply(seq(0, 70), function(i)
        probCounts2MN(300, 70, 110, i), numeric(1)))
    expect_equal(res, 1, tolerance=0.0001)
})

test_that("probCounts3N works", {
    res <- sum(vapply(seq(0, 20), function(k) probCounts3N(25, 20, 30, 140, k),
                      numeric(1)))
    expect_equal(res, 1, tolerance=0.0001)
})

