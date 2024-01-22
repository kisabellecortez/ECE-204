#include <iostream>
#include <functional>
#include <cmath>
#include <cassert>
#include "tlinalg.hpp"

// Function declarations
int main();
template <unsigned int n>
vec<n> markov_chain(
  matrix<n, n> A,
  vec<n> v0,
  double eps_step,
  unsigned int max_iterations
);

// Function definitions
int main() {
  ////////////////////////////////////////////
  // PROJECT
  // This is code that tests the project.
  ////////////////////////////////////////////

  vec<6> v0{ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  // You should have many more digits in the mantissa
  // than is shown here... Use 'format long' in Matlab
  matrix<6, 6> A{
    {0.3957, 0.1931, 0.0224, 0.8002, 0.4276, 1.0000},
    {0.8426, 0.4123, 0.9964, 0.3864, 0.6946, 0.0000},
    {0.7730, 0.7306, 0.1065, 0.3964, 0.9449, 0.0000},
    {0.2109, 0.7501, 0.4547, 0.7366, 0.3298, 0.0000},
    {0.6157, 0.8470, 0.4711, 0.3926, 0.8364, 0.0000},
    {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000}
  };

  // This should throw an exception
  try {
    std::cout << markov_chain<6>( A, v0, 1e-5, 100 )
              << std::endl;
  } catch ( std::invalid_argument &e ) {
    std::cout << "A is not stochastic" << std::endl;
  }

  // Make 'A' into a markov_chain matrix
  for ( unsigned int j{ 0 }; j < 6; ++j ) {
    double column_sum{ 0.0 };

    for ( unsigned int i{ 0 }; i < 6; ++i ) {
      column_sum += A( i, j );
    }

    for ( unsigned int i{ 0 }; i < 6; ++i ) {
      A( i, j ) /= column_sum;
    }
  }

  // This should print
  //  [0.139434 0.065835 0.010921 0.295037 0.132249;
  //   0.296910 0.140568 0.485788 0.142467 0.214827;
  //   0.272385 0.249088 0.051923 0.146154 0.292240;
  //   0.074316 0.255736 0.221686 0.271588 0.102001;
  //   0.216956 0.288773 0.229682 0.144753 0.258683]
  std::cout << A << std::endl;

  // This should print
  //     [0.123697 0.247392 0.202221 0.193653 0.233038]'
  std::cout << markov_chain<6>( A, v0, 1e-5, 100 )
            << std::endl;

  // Change 'A' so that the column sums are still 1.0,
  // but there is a negative entry in (0, 0).
  //  - Ethan Maeda noted that the second should be A( 1, 0 )
  A( 0, 0 ) -= 1.1;
  A( 1, 0 ) += 1.1;

  // This should throw an exception
  try {
    std::cout << markov_chain<6>( A, v0, 1e-5, 100 )
              << std::endl;
  } catch ( std::invalid_argument &e ) {
    std::cout << "A is not stochastic" << std::endl;
  }

  // PROJECT Question 5
  //
  // matrix<6, 6> B{ {...}, {...}, {...}, {...}, {...}, {...} };
  // Stochastic matrix
  // vec<6> u0{ 0.1, 0.3, 0.2, 0.1, 0.2, 0.1 };
  // Stochastic vector
  // std::cout << markov_chain<6>( B, u0, 1e-5, 1000 )
  //           << std::endl;

  matrix<6, 6> Q5{
    {4.9602e-02, 3.9491e-01, 8.0424e-01, 4.3817e-02, 1.3315e-01, 7.6256e-01},
    {9.8110e-01, 3.2924e-01, 9.5337e-01, 5.1886e-01, 6.0196e-01, 9.1649e-01},
    {6.6801e-01, 7.2262e-01, 6.1835e-01, 4.2310e-01, 1.2395e-01, 5.8291e-01},
    {8.7204e-01, 9.1476e-01, 1.7868e-01, 9.6011e-01, 4.2699e-01, 7.7038e-01},
    {7.4930e-01, 1.2133e-01, 8.5209e-01, 5.7775e-01, 9.3918e-01, 7.8593e-01},
    {2.3446e-03, 1.6840e-01, 6.4882e-01, 4.6606e-01, 9.4814e-01, 2.1560e-01}
  };

  return 0;
}

////////////////////////////////////////////
// PROJECT
// This is the function you need to
// implement
////////////////////////////////////////////

template <unsigned int n>
vec<n> markov_chain(
  matrix<n, n> A,
  vec<n> v0,
  double eps_step,
  unsigned int max_iterations
) {
  // Ensure that 'A' represents a stocastic matrix
  //  - All entries are non-negative
  //  - All of the rows add up to '1.0' with an
  //    allowed error of eps_step

    //check if stochastic matrix
    for(int k = 1; k <= n; ++k){ //loop through columns
      double sum = 0.0; //create variable to store sums of columns

      for(int j = 1; j <= n; ++j){ //loop through rows

        if(A(j, k) < 0){ //check if entry is negative
          throw std::invalid_argument("There is a negative number.")
        }

        sum += A(j, k)

      }

      if(std::abs(1.0 - sum) > n * eps_steps){ //check if column doesn't equal to 1
        throw std::invalid_argument("A column does not equal to 1.")
      }
    }

    

  // Iterate as necessary
  return vec<n>{};
}