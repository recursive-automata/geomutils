#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector clip(NumericVector x, double a, double b){
    return clamp(a, x, b);
}
