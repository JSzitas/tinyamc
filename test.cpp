//
// Created by jsco on 12/22/24.
//

//#include "heston.h"
#include "tinyamc.h"
#include "utils.h"


template <typename T>
void print_vector(const std::vector<T>& vec) {
  std::cout << "[";
  if (!vec.empty()) {
    std::for_each(vec.begin(), vec.end() - 1, [](const T& element) {
      std::cout << element << ", ";
    });
    std::cout << vec.back(); // Print the last element without a trailing comma
  }
  std::cout << "]" << std::endl;
}

template <typename T>
void print_square_matrix(const std::vector<T>& mat, size_t n) {
  if (mat.size() != n * n) {
    throw std::runtime_error("Matrix size is not consistent with dimensions.");
  }

  for (size_t i = 0; i < n; ++i) {
    std::cout << "[";
    for (size_t j = 0; j < n; ++j) {
      std::cout << mat[j * n + i] << (j < n - 1 ? ", " : ""); // Column-major access
    }
    std::cout << "]" << std::endl;
  }
  std::cout << std::endl; // Add extra newline for clarity
}

int main() {
  std::vector<double> x = {0.846155384215181, -0.62561632570142, 1.07094982024803, 1.06519898819057,
    -0.146513230426962, -1.52211873602427, -1.46057248176692, 1.3016834667554,
    -0.218233921680137, -0.0344228983329794};
  std::vector<double> y = {-0.0638114572918477, 1.3805276003709, 1.09576670963823, -0.705286887422382,
    0.798331909551031, 0.937896718947201, 1.33330195307445, -0.0362470387381375,
    -1.08536463844442, 0.0808272813021583};

  auto res = polynomial_regression(x, y, 3);
  //print_square_matrix(res, std::sqrt(res.size()));
  print_vector(res);
/*

  std::mt19937 gen(42);
  std::normal_distribution<> d{0, 1};

  std::vector<int> n_obs = {100, 200, 500,
                            1000, 2000, 5000,
                            10000, 20000, 50000,
                            100000}; // Reduced n_obs for testing
  int max_degree = 10;
  for (int n : n_obs) {
    std::cout << "n = " << n << std::endl;
    for (int degree = 0; degree <= max_degree; ++degree) {
      //std::cout << "\tDegree = " << degree << ": ";
      auto run_regression = [&]() {
        std::vector<double> x(n), y(n), coef(degree+1);
        for (int i = 0; i < n; ++i) {
          x[i] = d(gen);
          y[i] = d(gen); // Initialize y with random data
        }
        for(size_t i = 0; i < degree+1; ++i) {
          coef[i] = d(gen);
        }
        std::transform(x.begin(), x.end(), y.begin(), [&](double val) {
          double res = 0;
          double current_power = 1;
          for (int k = 0; k <= degree; k++){
            res += coef[k] * current_power;
            current_power *= val;
          }
          return res + d(gen);
        });
        std::vector<double> coefficients = polynomial_regression(x, y, degree);
        // test coefficients match
        bool fail = false;
        for(size_t i = 0; i < degree+1; ++i) {
          if(std::abs(coef[i] - coefficients[i]) > 0.05){
            if(!fail) {
              std::cout << "\tDegree = "<< degree << ": " << "Coefficients do not match" << std::endl;
            }
            std::cout << "coef: " << coef[i] << " | " << coefficients[i] << std::endl;
            fail = true;
          }
        }
        //if(fail) std::throw_with_nested(std::runtime_error("Coefficients do not match"));
      };
      run_regression();
    }
  }
  */
  return 0;
}


