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

  matrix<6, 6> B{
    {0.233393157265680, 0.071852318385612, 0.257150435427443, 0.187448140069635, 0.199716134899723, 0.356307680696386},
    {0.259481396638629, 0.141094995856625, 0.130399988888121, 0.227030303845137, 0.222963152938832, 0.016064514562381},
    {0.036377798363843, 0.247035999999946, 0.215001647744527, 0.155158087836934, 0.218664884513753, 0.139749781462134},
    {0.261653955070667, 0.248940472652941, 0.038118881796087, 0.008449918951417, 0.115411827108052, 0.023300491836811},
    {0.181151381097160, 0.040664049387831, 0.113309488666397, 0.200916733694215, 0.192872742485127, 0.049017762836316},
    { 0.027942311564020, 0.250412163717045, 0.246019557477425, 0.220996815602662, 0.050371258054513, 0.415559768605972}
  };

  vec<6> u0{ 0.1, 0.3, 0.2, 0.1, 0.2, 0.1 }; 

  std::cout << markov_chain<6>(B, u0, 1e-5, 1000)
            <<std::endl; 

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
          throw std::invalid_argument("There is a negative number.");
        };

        sum += A(j, k);

      };

      if(std::abs(1.0 - sum) > n * eps_step){ //check if column doesn't equal to 1
        throw std::invalid_argument("A column does not equal to 1.");
      };
    }

    for(int i = 1; i <= max_iterations; ++i){
      vec<n> vk = A * v0; 

      if((vk - v0).norm() < eps_step){
        return vk;
      };

      v0 = vk; 
    }

    throw std::runtime_error("Did not converge.");

  // Iterate as necessary
  return vec<n>{};
}