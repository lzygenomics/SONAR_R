#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

static double sonar_negloglik_impl(const NumericVector& x,
                                  const NumericMatrix& u,
                                  const NumericVector& N,
                                  const NumericMatrix& w,
                                  const NumericMatrix& y,
                                  int spot) {
  const int xindex = u.nrow();
  const int genes = u.ncol();
  const int spots = N.size();
  const double theta = x[xindex];

  if (!R_finite(theta) || theta <= 0.0) {
    return R_PosInf;
  }

  double result = 0.0;

  for (int n = 0; n < spots; ++n) {
    const double weight = w(spot, n);
    if (weight == 0.0) {
      continue;
    }

    double neighbor_loglik = 0.0;

    for (int j = 0; j < genes; ++j) {
      double nsum = 0.0;
      for (int k = 0; k < xindex; ++k) {
        nsum += x[k] * u(k, j);
      }

      const double mu = N[n] * nsum;
      const double count = y(n, j);

      if (!R_finite(nsum) || nsum <= 0.0 || !R_finite(mu) || mu <= 0.0) {
        return R_PosInf;
      }

      double first = 0.0;
      if (count > 0.0) {
        first = R::lgammafn(count + theta) - R::lgammafn(theta);
      }

      neighbor_loglik += first +
        theta * (std::log(theta) - std::log(theta + mu)) +
        count * (std::log(mu) - std::log(mu + theta));
    }

    result += neighbor_loglik * weight;
  }

  return -result;
}

// [[Rcpp::export]]
double sonar_negloglik_cpp(NumericVector x,
                           NumericMatrix u,
                           NumericVector N,
                           NumericMatrix w,
                           NumericMatrix y,
                           int spot) {
  return sonar_negloglik_impl(x, u, N, w, y, spot - 1);
}

// [[Rcpp::export]]
NumericVector sonar_grad_cpp(NumericVector x,
                             NumericMatrix u,
                             NumericVector N,
                             NumericMatrix w,
                             NumericMatrix y,
                             int spot) {
  const int xindex = u.nrow();
  const int genes = u.ncol();
  const int spots = N.size();
  const int spot0 = spot - 1;
  const double theta = x[xindex];

  NumericVector gradient(xindex + 1);

  if (!R_finite(theta) || theta <= 0.0) {
    std::fill(gradient.begin(), gradient.end(), R_NaN);
    return gradient;
  }

  double theta_result = 0.0;
  NumericVector beta_result(xindex);

  for (int n = 0; n < spots; ++n) {
    const double weight = w(spot0, n);
    if (weight == 0.0) {
      continue;
    }

    double theta_neighbor = 0.0;
    NumericVector beta_neighbor(xindex);

    for (int j = 0; j < genes; ++j) {
      double nsum = 0.0;
      for (int k = 0; k < xindex; ++k) {
        nsum += x[k] * u(k, j);
      }

      const double mu = N[n] * nsum;
      const double denom = theta + mu;
      const double count = y(n, j);

      if (!R_finite(nsum) || nsum <= 0.0 || !R_finite(mu) || mu <= 0.0 ||
          !R_finite(denom) || denom <= 0.0) {
        std::fill(gradient.begin(), gradient.end(), R_NaN);
        return gradient;
      }

      double first_theta = 0.0;
      if (count > 0.0) {
        first_theta = R::digamma(count + theta) - R::digamma(theta);
      }

      theta_neighbor += first_theta + std::log(theta) - std::log(denom) +
        1.0 - (theta + count) / denom;

      for (int k = 0; k < xindex; ++k) {
        const double ukj = u(k, j);
        beta_neighbor[k] +=
          -theta * N[n] * ukj / denom +
          count * (ukj / nsum - N[n] * ukj / denom);
      }
    }

    theta_result += theta_neighbor * weight;
    for (int k = 0; k < xindex; ++k) {
      beta_result[k] += beta_neighbor[k] * weight;
    }
  }

  for (int k = 0; k < xindex; ++k) {
    gradient[k] = -beta_result[k];
  }
  gradient[xindex] = -theta_result;

  return gradient;
}
