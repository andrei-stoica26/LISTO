test_that("pOverlapMNk works", {
    res <- sum(vapply(seq(0, 70), function(i)
        pOverlapMNk(300, 70, 110, i), numeric(1)))
    expect_equal(res, 1, tolerance=0.0001)
})

test_that("gene enrichment functions work", {
    m <- genesER(c('AURKA', 'TOP2A', 'CENPF', 'PTTG2', 'MKI67', 'BIRC5',
                   'RRM2'),
                 'human')
    expect_true('chromosome segregation' %in% m@result$Description)
    expect_equal(termGenes(m, 'chromosome segregation', 'meiosis I'),
                 c('BIRC5', 'CENPF', 'MKI67'))
})
