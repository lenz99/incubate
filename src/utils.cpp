#include <cpp11.hpp>
#include <Rmath.h>

using namespace cpp11;

//' Difference in Log-Space
//'
//' log-space difference of two values
//' @param lx logarithm of first value
//' @param ly logarithm of second value
//' @return log-space difference: log(exp(lx)-exp(ly))
//' @export
[[cpp11::register]]
double logspace_sub_cpp(double lx, double ly) {
  // use logspace_sub from Rmath
  return logspace_sub(lx, ly);
}

//' Difference in log-space
//'
//' log-space difference of first two values in a vector
//' @param lxy vector with two values, calculates 2nd - 1st value
//' @return log-space difference
//' @export
[[cpp11::register]]
double logspace_sub2_cpp(doubles lxy) {
  // use logspace_sub from Rmath
  return logspace_sub(lxy[1], lxy[0]);
}
