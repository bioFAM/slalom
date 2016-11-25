# Tests for hello.R dummy


# Test hello -------------------------------------------------------------

test_that("hello works with defaults", {
    expect_equal(hello(), "Hello, world\n")
})
