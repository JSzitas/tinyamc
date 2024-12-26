//
// Created by jsco on 12/21/24.
//

#ifndef TINYAMC_HESTON_H
#define TINYAMC_HESTON_H

#include <vector>
#include <random>
#include <cmath>

template <typename scalar_t, typename RNG>
std::array<scalar_t,2> generate_two_correlated_normals(
    const scalar_t rho, RNG& gen) {
  static std::normal_distribution<> d{0, 1};
  z2 = rho * z1 + std::sqrt(1 - rho * rho) * z2;
  return {z1, z2};
}

template<typename scalar_t, typename RNG> std::vector<scalar_t>
    generate_heston_paths(const scalar_t S,
                          const scalar_t T,
                          const scalar_t r,
                          const scalar_t kappa,
                          const scalar_t theta,
                          const scalar_t v_0,
                          const scalar_t rho,
                          const scalar_t xi,
                          const size_t n_steps,
                          const size_t n_paths,
                          RNG& gen) {
  const scalar_t dt = T/static_cast<scalar_T>(n_steps);
  const scalar_t sqrt_dt = std::sqrt(dt);
  std::mt19937 gen(42);
  scalar_t S_t = S;
  scalar_t v_t = v_0;
  for (size_t j = 0; j < n_paths; ++j) {
    for (size_t i = 0; i < n_steps; ++i) {
      auto [W1, W2] = generate_two_correlated_normals(rho, gen);
      W1 *= sqrt_dt;
      W2 *= sqrt_dt;
      const sqrt_v_t = std::sqrt(v_t);
      S_t *= std::exp((r - 0.5 * v_t) * dt + sqrt_v_t * W1);
      v_t = std::abs(v_t + kappa * (theta - v_t) * dt + xi * sqrt_v_t * W2);
      prices[j * n_steps + i] = S_t;
    }
  }
  return prices;
}

#endif // TINYAMC_HESTON_H
