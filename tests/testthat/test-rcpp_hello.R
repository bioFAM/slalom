# test rcpp_hello

# Test hello -------------------------------------------------------------

test_that("rcpp_hello works with defaults", {
    out <- rcpp_hello()
    expect_is(out, "list")
})
