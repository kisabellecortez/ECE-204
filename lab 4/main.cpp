#include <iostream>
#include <functional>
#include <cmath>
#include <cassert>
#include <tuple>
#include "tlinalg.hpp"

// Function declarations
int main();

std::tuple<double *, vec<3> *, vec<3> *> rk4(
    std::function<vec<3>(double t, vec<3> y)> f,
    std::pair<double, double> t_rng, vec<3> y0,
    unsigned int n);

vec<3> ivp_evaluate(
    double t,
    std::tuple<double *, vec<3> *, vec<3> *> y,
    unsigned int n);

vec<3> ivp_spline_2pt(
    double t,
    double ts[2],
    vec<3> ys[2],
    vec<3> dys[2]);

void print(
    std::tuple<double *, vec<3> *, vec<3> *> out,
    unsigned int n);

vec<3> f1(double t, vec<3> y);
vec<3> f2(double t, vec<3> y);
vec<3> f3(double t, vec<3> y);

// Function definitions
int main()
{

    std::cout << f1(0.3, {1.0, 2.0, 3.0}) << std::endl;

    // Showing your code works
    std::cout.precision(8);

    auto result{rk4(
        f2, std::make_pair(0.0, 4.0), {1.0, 2.0, 3.0}, 8)};

    print(result, 8);

    delete[] std::get<0>(result);
    delete[] std::get<1>(result);
    delete[] std::get<2>(result);

    // Showing Runge Kutta is O(h^4)
    std::cout.precision(16);
    // Replace this with the solution found
    // in Matlab for y(10)
    vec<3> soln{-0.047714, -0.014473, 0.015358};

    for (unsigned int n{1}; n <= 128; n *= 2)
    {
        auto result{rk4(
            f3, std::make_pair(0.0, 10.0), {1.0, 2.0, 3.0}, n)};

        std::cout << norm(std::get<1>(result)[n] - soln)
                  << std::endl;

        delete[] std::get<0>(result);
        delete[] std::get<1>(result);
        delete[] std::get<2>(result);
    }

    // Showing your approximations using cubic splines
    std::cout.precision(8);

    auto result2a{rk4(
        f3, std::make_pair(0.0, 1.0), {1.0, 2.0, 3.0}, 8)};

    print(result2a, 8);

    auto result2b{rk4(
        f3, std::make_pair(0.0, 1.0), {1.0, 2.0, 3.0}, 4)};

    print(result2b, 4);

    for (unsigned int k{0}; k <= 16; ++k)
    {
        std::cout << "t = " << 0.0625 * k << std::endl;
        std::cout << "\t\ty = "
                  << ivp_evaluate(0.0625 * k, result2a, 8)
                  << std::endl;
        std::cout << "\t\ty = "
                  << ivp_evaluate(0.0625 * k, result2b, 4)
                  << std::endl;
    }

    delete[] std::get<0>(result2a);
    delete[] std::get<1>(result2a);
    delete[] std::get<2>(result2a);

    delete[] std::get<0>(result2b);
    delete[] std::get<1>(result2b);
    delete[] std::get<2>(result2b);

    return 0;
}

// The systems of differential equations
vec<3> f1(double t, vec<3> y)
{
    return vec<3>{
        y(0) + y(1) - 3 * y(2) + 1,
        2 * y(0) + y(1) - y(2) + 2,
        y(0) + 4 * y(1) - y(2) - 1};
}

vec<3> f2(double t, vec<3> y)
{
    return vec<3>{
        -0.3 * y(0) + 0.09 * y(1) - 0.11 * y(2),
        -0.7 * y(1) - 0.08 * y(0) + 0.03 * y(2),
        -0.4 * y(2) + 0.02 * y(0) - 0.15 * y(1)};
}

vec<3> f3(double t, vec<3> y)
{
    return vec<3>{
        -0.49 * y(0) + 0.89 * y(1) - 0.03 * y(2),
        -0.79 * y(0) - 0.43 * y(1) - 0.45 * y(2),
        0.03 * y(0) + 0.56 * y(1) - 0.38 * y(2)};
}

// These are the non-linear ordinary
// differential equations that define
// the Lorenz equations. You can
// give these any initial condition you
// wish, but ultiplate, if the outputs
// are plotted as coordinates in 3d,
// they form the butterfly found here:
//   https://en.wikipedia.org/wiki/Lorenz_system
vec<3> lorenz(double t, vec<3> y)
{
    double rho{28.0};
    double sigma{10.0};
    double beta{8.0 / 3.0};

    return vec<3>{
        sigma * (y(1) - y(0)),
        y(0) * (rho - y(2)) - y(1),
        y(0) * y(1) - beta * y(2)};
}

////////////////////////////////////////////
// PROJECT
// This is the function you need to
// implement
////////////////////////////////////////////

std::tuple<double *, vec<3> *, vec<3> *> rk4(
    std::function<vec<3>(double t, vec<3> y)> f,
    std::pair<double, double> t_rng, vec<3> y0,
    unsigned int n)
{
    /*
    // This is the source code for the RK4 method
    // for a single initial-value problem.
    assert( n > 0 );

    //                              tf - t0
    // Calculate the step size h = ---------
    //                                 n
    double h{ (t_rng.second - t_rng.first)/n };

    // Allocate memory for the three output arrays
    double *ts{ new double[n + 1] };
    double *ys{ new double[n + 1] };
    double *dys{ new double[n + 1] };

    // Assign the initial condtions
    ts[0] = t_rng.first;
    ys[0] = y0;
    dys[0] = f( ts[0], ys[0] );

    for ( unsigned int k{0}; k < n; ++k ) {
      ts[k + 1] = ts[0] + h*(k + 1);

      // Calculate the four slopes
      double s0{ dys[k] };
      double s1{ f( ts[k] + h/2.0, ys[k] + h*s0/2.0 ) };
      double s2{ f( ts[k] + h/2.0, ys[k] + h*s1/2.0 ) };
      double s3{ f( ts[k + 1], ys[k] + h*s2 ) };

      // Find the approximation of y(t )
      //                              k
      ys[k + 1] = ys[k] + h*(
        s0 + 2.0*s1 + 2.0*s2 + s3
      )/6.0;
      // Store the slope y'(t ) = f( t , y  )
      //                     k        k   k
      dys[k + 1] = f( ts[k + 1], ys[k + 1] );
    }
    return std::make_tuple( ts, ys, dys );
    */
     double h{ (t_rng.second - t_rng.first)/n };

     // Allocate memory for the three output arrays
     double *ts{ new double[n + 1] };
     vec<3> *ys{ new vec<3>[n + 1] };
      vec<3> *dys{ new vec<3>[n + 1] };

     // Assign the initial condtions
     ts[0] = t_rng.first;
     ys[0] = y0;
     dys[0] = f( ts[0], ys[0] );

     for ( unsigned int k{0}; k < n; ++k ) {
       ts[k + 1] = ts[0] + h*(k + 1);

       // Calculate the four slopes
       vec<3> s0{ dys[k] };
       vec<3> s1{ f( ts[k] + h/2.0, ys[k] + h*s0/2.0 ) };
       vec<3> s2{ f( ts[k] + h/2.0, ys[k] + h*s1/2.0 ) };
       vec<3> s3{ f( ts[k + 1], ys[k] + h*s2 ) };

       // Find the approximation of y(t )
       //                              k
       ys[k + 1] = ys[k] + h*(
         s0 + 2.0*s1 + 2.0*s2 + s3
       )/6.0;
       // Store the slope y'(t ) = f( t , y  )
       //                     k        k   k
       dys[k + 1] = f( ts[k + 1], ys[k + 1] );
     }

      return std::make_tuple(ts, ys, dys);
}

vec<3> ivp_evaluate(
    double t,
    std::tuple<double *, vec<3> *, vec<3> *> y,
    unsigned int n)
{
    double *ts = std::get<0>(y);
    vec<3> *ys = std::get<1>(y);
    vec<3> *dys = std::get<2>(y);

    bool isFound = false;
    unsigned int k = 0;

    // Find the index k such that ts[k-1] < t < ts[k]
    for (unsigned int i = 1; i <= n; ++i)
    {
        if (t >= ts[i - 1] && t <= ts[i])
        {
            k = i;
            isFound = true;
            break;
        }
    }

    if (!isFound)
    {
        throw std::domain_error("t-value is outside the interval [t0, tf]");
    }

    // Check if t-value is one of the n+1 entries of the array ts
    if (t == ts[k - 1])
    {
        return ys[k - 1];
    }
    else
    {
        // Evaluate cubic splines at the point t
        double t_values[2] = {ts[k - 1], ts[k]};
        vec<3> y_values[2] = {ys[k - 1], ys[k]};
        vec<3> dy_values[2] = {dys[k - 1], dys[k]};

        return ivp_spline_2pt(t, t_values, y_values, dy_values);
    }
}

void print(
    std::tuple<double *, vec<3> *, vec<3> *> out,
    unsigned int n)
{
    // Print the n + 1 t-values
    std::cout << "The " << (n + 1) << " t-values:" << std::endl;

    std::cout << std::get<0>(out)[0];

    for (std::size_t k{1}; k <= n; ++k)
    {
        std::cout << ", " << std::get<0>(out)[k];
    }

    std::cout << std::endl
              << std::endl;
    std::cout << "The approximation of the three functions:"
              << std::endl;

    // For each of the three functions, print the corresponding
    // approximations for the above n + 1 t-values.
    for (std::size_t i{0}; i < 3; ++i)
    {
        std::cout << std::get<1>(out)[0](i);

        for (std::size_t k{1}; k <= n; ++k)
        {
            std::cout << ", " << std::get<1>(out)[k](i);
        }

        std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << "The slopes at the approximations:"
              << std::endl;

    for (std::size_t i{0}; i < 3; ++i)
    {
        std::cout << std::get<2>(out)[0](i);

        for (std::size_t k{1}; k <= n; ++k)
        {
            std::cout << ", " << std::get<2>(out)[k](i);
        }

        std::cout << std::endl;
    }

    std::cout << std::endl;
}

vec<3> ivp_spline_2pt(
    double t,
    double ts[2],
    vec<3> ys[2],
    vec<3> dys[2])
{
    double h = ts[1] - ts[0];
    double delta = (t - ts[0]) / h;
    assert((0.0 <= delta) && (delta <= 1.0));

    return (((h * (dys[0] + dys[1]) + 2.0 * (ys[0] - ys[1])) * delta - (h * (2.0 * dys[0] + dys[1]) + 3.0 * (ys[0] - ys[1]))) * delta + h * dys[0]) * delta + ys[0];
}