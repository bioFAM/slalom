# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#' Hello!
#'
#' A function that says hello.
#'
#' @param txt character string to be added to "Hello, " string.
#'
#' @export
#' @examples
#' hello("world")
hello <- function(txt = "world") {
    out <- paste0("Hello, ", txt, "\n")
    cat(out)
    out
}

