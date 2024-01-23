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
    {3.084625352677238e-01, 5.479347904723656e-01, 1.028229575076087e-01, 2.362337358920694e-01, 4.237720030896696e-01, 6.547064913012779e-01},
    {8.516480280975659e-01, 6.280441694128458e-01, 1.395781838033429e-01, 2.394217037762256e-01, 4.438425244122398e-01, 5.167364085231622e-01},
    {9.393170825177386e-01, 1.927750992433191e-01, 3.152049894797712e-01, 4.250708761523204e-01, 1.224363687105852e-01, 6.953524952537744e-01},
    {6.886279494317848e-01, 8.141847382630821e-01, 7.532793962721368e-01, 3.248275067488126e-01, 2.559292387453460e-02, 2.931359721403338e-01},
    {6.623268198900989e-02, 1.590382026051161e-01, 4.822004988781755e-01, 2.188593696826540e-01, 4.065125736025650e-01, 5.934552929353537e-01},
    {2.251602330132394e-01, 9.970123770177969e-01, 1.780947165125635e-02, 3.695431441239660e-01, 2.014672781172335e-01, 9.463154017778177e-01}
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
          throw std::invalid_argument("There is a negative number.")
        }

        sum += A(j, k)

      }

      if(std::abs(1.0 - sum) > n * eps_steps){ //check if column doesn't equal to 1
        throw std::invalid_argument("A column does not equal to 1.")
      }
    }

    for(int i = 1; i <= max_iterations; ++i){
      vec<n> vk = A * v0; 

      if((vk - v0).norm() < eps_step){
        return vk
      }

      v0 = vk; 
    }

    throw std::runtime_error("Did not converge.")

  // Iterate as necessary
  return vec<n>{};
}