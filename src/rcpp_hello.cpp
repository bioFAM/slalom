#include <Rcpp.h>
using namespace Rcpp;

// This is a simple function using Rcpp that creates an R list
// containing a character vector and a numeric vector.
//
// Learn more about how to use Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//
// and browse examples of code using Rcpp at:
//
//   http://gallery.rcpp.org/
//

//' @title
//' rcpp_hello
//' @description
//' This is a simple function using Rcpp that creates an R list containing a
//' character vector and a numeric vector.
//'
//' @details
//' Not much to report.
//'
//' @return
//' a list with a character vector and a numeric vector
//'
//' @useDynLib slalom
//' @importFrom Rcpp evalCpp
//' @export
//' @examples
//' rcpp_hello()
// [[Rcpp::export]]
List rcpp_hello() {
  CharacterVector x = CharacterVector::create("foo", "bar");
  NumericVector y   = NumericVector::create(0.0, 1.0);
  List z            = List::create(x, y);
  return z;
}
