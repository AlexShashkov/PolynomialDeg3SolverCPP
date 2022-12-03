#ifndef POLYROOTS1234_HPP
#define POLYROOTS1234_HPP
#define PR_NUMBERS_OF_ROOTS_EQUAL     0
#define PR_AT_LEAST_ONE_ROOT_LOST    -1
#define PR_AT_LEAST_ONE_ROOT_IS_FAKE -2
#define PR_2_INFINITE_ROOTS          -3

#define PR_DISCRIMINANT_USE_TRADITIONAL_OPERATIONS_NORMALIZED 0
#define PR_DISCRIMINANT_USE_TRADITIONAL_OPERATIONS_NORMALIZED 0
#define PR_DISCRIMINANT_USE_FMA_UNNORMALIZED 0
#define PR_DISCRIMINANT_USE_FMA_NORMALIZED 0
//#include <array>
#include <cmath>
#include <numbers> // std::numbers::pi_v<fp_t>, requires -std=c++20
#include <chrono> // for testing only
#include <random> // for testing only
#include <cstdlib> // for testing only
#include <cassert>
#include <iostream>
#include <iomanip>
#include <complex>

// checks attainable number of real roots in a polynomial: a*x^4 + b*x^3 + c*x^2 + d*x + e; multiple root is treated as separate roots
template<typename fp_t>
int number_of_roots(unsigned P, // polynomial degree
                    fp_t a, fp_t b, fp_t c, fp_t d, fp_t e); // polynomial coefficients
// Compares two vectors of roots; root orderings play no role. For each entry in (roots_ground_truth),
// the closest entry in (roots_to_check) is found and corresponding distance found. Among such distances
// the largest will be stored to (max_deviation)
template<typename fp_t>
int compare_roots(
        unsigned N_roots_to_check, // number of roots in (roots_to_check)
        unsigned N_roots_ground_truth,  // number of roots in (roots_ground_truth)
        std::vector<fp_t> &roots_to_check, // one should take into account only first (N_roots_to_check) roots here
        std::vector<fp_t> &roots_ground_truth, // one should take into account only first (N_roots_ground_truth) roots here
        fp_t &max_absolute_error, // here the greatest among the smallest deviations of the roots in (roots_to_check) and (roots_ground_truth)
        // will be placed
        // here the greatest relative error among all the roots found will be placed
        fp_t &max_relative_error);


// Creates a test polynomial, both in the form of roots, e.g. (x-roots[0])*(x-roots[1])*(quadratic polynomial with no real roots) and
// represented by coefficients, e.g. (coefficients[4]=1)*x^4 + coefficients[3]*x^3 + coefficients[2]*x^2 + coefficients[1]*x + coefficients[0].
// The highest-degree coefficient always equals 1. The function returns the actual number of different real roots placed into the vector
// (roots) (complex roots are not placed there). Negative return values may mean internal implementation error
template<typename fp_t>
int generate_polynomial(
        unsigned P, // polynomial degree
        unsigned N_pairs_of_complex_roots, // how many pairs of complex conjugate roots to introduce
        unsigned N_clustered_roots, // how many clustered roots to introduce; all the clustered roots are real
        unsigned N_multiple_roots, // how many multiple roots to introduce; all multiple roots are real
        fp_t max_distance_between_clustered_roots, // maximal distance between the closest of the clustered roots
        fp_t root_sweep_low,
        fp_t root_sweep_high, // low and high boundaries of real roots; imaginary parts of complex conjugate roots are in the same range
        std::vector<fp_t> &roots, // storage where to put the roots; size should exceed P-1
        std::vector<fp_t> &coefficients);

template<typename fp_t>
int compare_roots_complex(unsigned N_roots_to_check, // number of roots in roots_to_check
                          unsigned N_roots_ground_truth,  // number of roots in roots_ground_truth
                          std::vector<std::complex<fp_t>> &roots_to_check, // one should take into account only first (N_roots_to_check) rots
                          std::vector<fp_t> &roots_ground_truth, // one should take into account only first (N_roots_ground_truth) rots
                          fp_t &max_absolute_error, // here the greatest among the smallest deviations of the roots in (roots_to_check) and (roots_ground_truth)
        // will be placed
        // here the greatest relative error among all the roots found will be placed
                          fp_t &max_relative_error);

#endif //POLYROOTS1234_HPP
