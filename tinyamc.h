
#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>

#include "tinyqr.h"

using date = double;

#include <algorithm>
#include <armadillo>
#include <cmath>
#include <stdexcept>
#include <vector>

template <typename scalar_t>
void transpose(std::vector<scalar_t>& A, const size_t n) {
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = i + 1; j < n; ++j) {
      std::swap(A[i * n + j], A[j * n + i]);
    }
  }
}

template <typename scalar_t>
void transpose(std::vector<scalar_t>& A) {
  const size_t n = std::sqrt(A.size());
  transpose(A, n);
}

template <typename scalar_t>
void invert_symmetric_inplace(std::vector<scalar_t>& A, const size_t n) {
  // Perform Cholesky decomposition: A = L * L^T
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j <= i; ++j) {
      scalar_t sum = A[i * n + j];
      for (size_t k = 0; k < j; ++k) {
        sum -= A[i * n + k] * A[j * n + k];
      }
      if (i == j) {
        A[i * n + j] = std::sqrt(sum);
      } else {
        A[i * n + j] = sum / A[j * n + j];
      }
    }
  }
  for (size_t i = 0; i < n; ++i) {
    A[i * n + i] = 1.0 / A[i * n + i];
    for (size_t j = i + 1; j < n; ++j) {
      scalar_t sum = 0;
      for (size_t k = i; k < j; ++k) {
        sum -= A[j * n + k] * A[k * n + i];
      }
      A[j * n + i] = sum / A[j * n + j];
    }
  }
  // Compute A_inverse = L^T_inverse * L_inverse in-place
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = i; j < n; ++j) {
      scalar_t sum = 0;
      for (size_t k = j; k < n; ++k) {
        sum += A[k * n + i] * A[k * n + j];
      }
      A[i * n + j] = sum;
    }
  }
}

template <typename scalar_t, const bool only_itm = false>
std::vector<scalar_t> polynomial_regression(
    const std::vector<scalar_t>& x, const std::vector<scalar_t>& y,
    const size_t degree) {
  const auto n = x.size();
  const auto n_coeffs = degree + 1;
  std::vector<scalar_t> XTX(n_coeffs * n_coeffs, 0), XTY(n_coeffs, 0);
  std::vector<scalar_t> power_sums(2 * degree + 1, 0);
  power_sums[0] = static_cast<scalar_t>(n);
  for (size_t i = 0; i < n; ++i) {
    // skip oom paths
    if constexpr (only_itm) {
      if(x[i] <= 0) {
        power_sums[0]--;
        continue;
      }
    }
    scalar_t x_pow_j = x[i];
    XTY[0] += y[i];
    for (size_t j = 1; j <= 2 * degree; ++j) {
      if (j < n_coeffs)
        XTY[j] += x_pow_j * y[i];
      power_sums[j] += x_pow_j;
      x_pow_j *= x[i];
    }
  }

  for (size_t j = 0; j < n_coeffs; ++j) {
    for (size_t k = 0; k <= j; ++k) {
      // write to both lower and upper diagonal in single pass
      XTX[j * n_coeffs + k] = power_sums[j + k];
      XTX[k * n_coeffs + j] = power_sums[j + k];
    }
  }
  //return XTX;
  invert_symmetric_inplace(XTX, n_coeffs);
  // compute coefficients beta
  power_sums.resize(n_coeffs);
  std::fill(power_sums.begin(), power_sums.end(), 0);
  for (size_t i = 0; i < n_coeffs; ++i) {
    for (size_t j = 0; j < n_coeffs; ++j) {
      // reuse as coefficients :)
      // if we access XTX in upper triangle, get lower triangle value instead
      if(i > j) {
        power_sums[i] += XTX[j * n_coeffs + i] * XTY[j];
        continue;
      }
      power_sums[i] += XTX[i * n_coeffs + j] * XTY[j];
    }
  }

  // compute hat vector
  std::vector<scalar_t> hat(n, 0);
  std::vector<scalar_t> Xi(n_coeffs, 1.0);
  std::vector<scalar_t> Xi_XTX(n_coeffs, 0.0);

  for (size_t i = 0; i < n; ++i) {
    // Construct the i-th row of X
    scalar_t x_pow = x[i];
    for (size_t j = 1; j < n_coeffs; ++j) { // Fill the rest of the Xi row as powers of x
      Xi[j] = x_pow;
      x_pow *= x[i];
    }
    // Compute Xi * XTX
    for (size_t j = 0; j < n_coeffs; ++j) {
      for (size_t k = 0; k < n_coeffs; ++k) {
        // adjustment because we only have lower diagonal values
        if(j > k) {
          Xi_XTX[j] += Xi[k] * XTX[k * n_coeffs + j];
          continue;
        }
        Xi_XTX[j] += Xi[k] * XTX[j * n_coeffs + k];
      }
    }

    // Compute (Xi * XTX) * (transpose(Xi)), which is just the dot product
    scalar_t hat_element = 0.0;
    for (size_t j = 0; j < n_coeffs; ++j) {
      hat_element += Xi_XTX[j] * Xi[j];
    }
    // Store the result in the hat vector
    hat[i] = hat_element;
    std::fill(Xi.begin()+1, Xi.end(), 1.0);
    std::fill(Xi_XTX.begin(), Xi_XTX.end(), 0.0);
  }
  // compute adjusted continuation values

  return hat;
}
/*
template <typename scalar_t>
std::vector<scalar_t> polynomial_regression_itm(
    const std::vector<scalar_t>& x, const std::vector<scalar_t>& y,
    const size_t degree) {
  const auto n = x.size();
  const auto n_coeffs = degree + 1;
  std::vector<scalar_t> XTX(n_coeffs * n_coeffs, 0), XTY(n_coeffs, 0);
  std::vector<scalar_t> power_sums(2 * degree + 1, 0);
  for (size_t i = 0; i < n; ++i) {
    // skip OOM paths
    if(x[i] <= 0) continue;
    scalar_t x_pow_j = 1;
    for (size_t j = 0; j <= 2*degree; ++j) {
      if (j < n_coeffs) XTY[j] += x_pow_j * y[i];
      power_sums[j] += x_pow_j;
      x_pow_j *= x[i];
    }
  }

  for (size_t j = 0; j < n_coeffs; ++j) {
    for (size_t k = 0; k <= j; ++k) {
      // write to both lower and upper diagonal in single pass
      XTX[j * n_coeffs + k] = power_sums[j + k];
      XTX[k * n_coeffs + j] = power_sums[j + k];
    }
  }
  return tinyqr::lm(XTX, XTY);
}
// returns value of continuation (i.e. fitted value) rather than coefficients
template <typename scalar_t>
std::vector<scalar_t> polynomial_regression_itm_fitted(
    const std::vector<scalar_t>& x, const std::vector<scalar_t>& y,
    const size_t degree) {
  const auto n = x.size();
  const auto n_coeffs = degree + 1;
  std::vector<scalar_t> XTX(n_coeffs * n_coeffs, 0), XTY(n_coeffs, 0);
  std::vector<scalar_t> power_sums(2 * degree + 1, 0);
  for (size_t i = 0; i < n; ++i) {
    // skip OOM paths
    if(x[i] <= 0) continue;
    scalar_t x_pow_j = 1;
    for (size_t j = 0; j <= 2*degree; ++j) {
      if (j < n_coeffs) XTY[j] += x_pow_j * y[i];
      power_sums[j] += x_pow_j;
      x_pow_j *= x[i];
    }
  }

  for (size_t j = 0; j < n_coeffs; ++j) {
    for (size_t k = 0; k <= j; ++k) {
      // write to both lower and upper diagonal in single pass
      XTX[j * n_coeffs + k] = power_sums[j + k];
      XTX[k * n_coeffs + j] = power_sums[j + k];
    }
  }
  const auto coef = tinyqr::lm(XTX, XTY);
  // compute hat matrix
  //X %*% coef %*% Xt;

}
*/

template<typename scalar_t> struct DiscounterFixed{
  scalar_t df(date i, date j) {
    // assumes fixed step - this is mainly to simplify implementation
    return 0.0995;
  }
};

template<typename scalar_t, typename MTM, typename Discounter> scalar_t amc(
    MTM& mtm,
    Discounter& disc,
    const size_t degree = 3,
    const size_t n_paths = 100,
    const size_t n_steps = 100) {
  // step 1 - generate MTM at each time step; this is actually something we
  // kind of ignore, since mtm is supposed to be providing those
  // we dont really care if they are pre-generated, or computed on the fly
  // we can thus safely proceed to the backward induction step; very fun
  // we assume this exists
  // use CV, PV for continuation value, present value
  const auto CV = mtm.last();
  for(size_t step = n_steps-1; step > 0; --step) {
    const auto PV = mtm(step);
    // carry out regression; note that we only want to select in the money paths
    polynomial_regression_itm(PV, CV, degree);



    // mtm should be of format


  }

}


