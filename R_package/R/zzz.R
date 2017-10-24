## For R 2.15.1 and later this also works. Note that calling loadModule() triggers
## a load action, so this does not have to be placed in .onLoad() or evalqOnLoad().
Rcpp::loadModule("SlalomModel", TRUE)

utils::globalVariables(c("n_gain", "term", "n_loss"), "slalom")

