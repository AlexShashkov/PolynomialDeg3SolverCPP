// g++ framework.cpp Methods.h -DMETHOD="Baydoun<fp_t> Solver" -DNTESTS=100000 -o Baydoun -std=c++20
// g++ framework.cpp Methods.h -DMETHOD="Vieta<fp_t> Solver" -DNTESTS=100000 -o Vieta -std=c++20

#ifndef POLYROOTS1234_HPP
#define POLYROOTS1234_HPP
#define PR_NUMBERS_OF_ROOTS_EQUAL     0
#define PR_AT_LEAST_ONE_ROOT_LOST    -1
#define PR_AT_LEAST_ONE_ROOT_IS_FAKE -2
#define PR_2_INFINITE_ROOTS          -3

#include <cmath>
#include <numbers> // std::numbers::pi_v<fp_t>, requires -std=c++20
#include <chrono> // for testing only
#include <random> // for testing only
#include <cstdlib> // for testing only
#include <cassert>
#include <iostream>
#include <iomanip>

#include "Methods.h"
using namespace implementations;

// algorithm selector in polyroots2()
#define PR_DISCRIMINANT_USE_TRADITIONAL_OPERATIONS_NORMALIZED 1
#define PR_DISCRIMINANT_USE_FMA_UNNORMALIZED                  2
#define PR_DISCRIMINANT_USE_FMA_NORMALIZED                    3

#ifndef METHOD
#define METHOD Baydoun<fp_t> Solver
#endif

#ifndef NTESTS
#define NTESTS 10000
#endif

int pr_discriminant_computation_method=0;
int pr_lost_roots_saving_1=0;
int pr_lost_roots_saving_2=0;
unsigned long long N_additive_feedback_fired=0ULL;
unsigned long long N_multiplicative_feedback_fired=0ULL;


/* computes (a*b - c*d) with precision not worse than 1.5*(unit of least precision) suggested in Claude-Pierre Jeannerod,
Nicolas Louvet, and Jean-Michel Muller, "Further Analysis of Kahan's Algorithm for the Accurate Computation of 2x2 Determinants".
Mathematics of Computation, Vol. 82, No. 284, Oct. 2013, pp. 2245-2264 */
template <typename fp_t> inline fp_t pr_product_difference(fp_t a, fp_t b, fp_t c, fp_t d)
{
  auto tmp = d * c;
  return fma(a, b, -tmp) + fma(-d, c, tmp);
}


// Creates a test polynomial, both in the form of roots, e.g. (x-roots[0])*(x-roots[1])*(quadratic polynomial with no real roots)
// as well as represented by its coefficients, e.g.
// (coefficients[4]=1)*x^4 + coefficients[3]*x^3 + coefficients[2]*x^2 + coefficients[1]*x + coefficients[0].
// The highest-degree coefficient always equals 1. The function returns the actual number of different real roots placed into the vector
// (roots) (complex roots are not placed there). Negative return values may mean internal implementation error
template<typename fp_t> int generate_polynomial(
unsigned P, // polynomial degree
unsigned N_pairs_of_complex_roots, // how many pairs of complex conjugate roots to introduce
unsigned N_clustered_roots, // how many clustered roots to introduce; all the clustered roots are real
unsigned N_multiple_roots, // how many multiple roots to introduce; all multiple roots are real
fp_t max_distance_between_clustered_roots, // maximal distance between the closest of the clustered roots
fp_t root_sweep_low, fp_t root_sweep_high, // low and high boundaries of real roots; imaginary parts of complex conjugate roots are in the same range
std::vector<fp_t> &roots, // storage where to put the roots; size should exceed P-1
std::vector<fp_t> &coefficients) // storage where to put the coefficients; size should exceed P
{
int n_simple_roots=P-2*N_pairs_of_complex_roots-N_clustered_roots-N_multiple_roots;
assert(N_clustered_roots!=1); assert(N_multiple_roots!=1); assert(n_simple_roots>=0); assert(P>0 && P<=4);
assert(max_distance_between_clustered_roots>static_cast<fp_t>(0.0L));
assert(root_sweep_high-root_sweep_low>2*P*max_distance_between_clustered_roots);

unsigned long long seed=std::chrono::system_clock::now().time_since_epoch().count()+std::rand(); // counts milliseconds
std::mt19937_64 rng(seed); // randomize seed from the clock
std::uniform_real_distribution<fp_t> rnr(root_sweep_low, root_sweep_high); // uniform random data generator for single roots
std::uniform_real_distribution<fp_t> rnc(static_cast<fp_t>(0.0L), max_distance_between_clustered_roots); // uniform random data generator for root clusters
fp_t re, im, u, v, root_mid_sweep=root_sweep_low+0.5*(root_sweep_high-root_sweep_low);
long double RE, IM, U, V, TMP; // high-precisioon counterparts of re, im, u, v

coefficients[P]=static_cast<fp_t>(1.0L); // invariant
switch (P)
  {
  case 0:
    coefficients[0]=rnr(rng); return 0;
  case 1:
    coefficients[0]=-(roots[0]=rnr(rng)); return 1;
  case 2:
    {
    if (N_pairs_of_complex_roots==1) // no real roots
      {
      re=rnr(rng); while ((im=rnr(rng))==static_cast<fp_t>(0.0L)) {}
      RE=re; IM=im;
      coefficients[1]=static_cast<fp_t>(-2.0L*RE); // -2*re
      coefficients[0]=static_cast<fp_t>(pr_product_difference(RE, RE, -IM, IM)); // re*re+im*im
      return 0;
      }
    else if (N_clustered_roots==2) // 2 close but distinct roots
      { roots[0]=re=rnr(rng); while ((im=rnc(rng))==static_cast<fp_t>(0.0L)) {} roots[1]=im=(re>=root_mid_sweep ? re-im : re+im); }
    else if (N_multiple_roots==2) // double root counted as a single root
      { roots[1]=roots[0]=im=re=rnr(rng); }
    else // 2 distinct single roots
      { roots[0]=re=rnr(rng); while ((im=rnr(rng))==re) {} roots[1]=im=rnr(rng); }
    RE=re; IM=im;

    coefficients[1]=static_cast<fp_t>( -RE-IM );
    coefficients[0]=static_cast<fp_t>( RE*IM );
    return 2; // return ((re!=im) ? 2 : 1);
    } // P=2
  case 3:
    {
    if (N_pairs_of_complex_roots==1) // one real root
      {
      re=rnr(rng); while ((im=rnr(rng))==static_cast<fp_t>(0.0L)) {} roots[0]=u=rnr(rng);
      RE=re; IM=im; U=u;
  
      IM=pr_product_difference(RE, RE, -IM, IM); // re*re+im*im
      RE*=-2.0L; // irreducible quadratic polynomial is (x^2 + re*x + im); multiply it by (x-u)
      coefficients[0]=static_cast<fp_t>(-IM*U); coefficients[2]=static_cast<fp_t>(RE-U);
      coefficients[1]=static_cast<fp_t>(std::fma(-RE,U,IM)); // im-re*u; 
      return 1;
      }
    else if (N_clustered_roots==3) // 3 clustered distinct roots
      {
      roots[0]=re=rnr(rng); while ((im=rnc(rng))==static_cast<fp_t>(0.0L)) {} while((u=rnc(rng))==static_cast<fp_t>(0.0L)) {}
      roots[1]=im=(re>root_mid_sweep ? roots[2]=u=(re-im-u), re-im : roots[2]=u=(re+im+u), re+im);
      }
    else if (N_clustered_roots==2) // 2 clustered roots, 1 single root; all distinct
      {
      roots[0]=re=rnr(rng); while((im=rnc(rng))==static_cast<fp_t>(0.0L)) {}
      roots[1]=im=(re>root_mid_sweep ? re-im : re+im); do { roots[2]=u=rnr(rng); } while (u==re || u==roots[1]);
      }
    else if (N_multiple_roots==3) // triple root counted as a single root
      { roots[2]=roots[1]=roots[0]=u=im=re=rnr(rng); }
    else if (N_multiple_roots==2) // double root and 1 single root; totally 2 roots
      { roots[1]=roots[0]=im=re=rnr(rng); while ((roots[2]=u=rnr(rng))==re) {} }
    else // 3 distinct single roots
      { roots[0]=re=rnr(rng); while ((roots[1]=im=rnr(rng))==re) {} do { roots[2]=u=rnr(rng); } while(u==re || u==im); }
    RE=re; IM=im; U=u;
    coefficients[2]=static_cast<fp_t>(-RE-IM-U); coefficients[0]=static_cast<fp_t>(-RE*IM*U);
    V=pr_product_difference(RE,IM,-RE,U); coefficients[1]=static_cast<fp_t>(std::fma(IM,U,V)); // re*im+re*u+im*u=im*u+(re*im-(-re*u));
    // if (re!=im && im!=u && u!=re) return 3;
    // if (re==im && im==u) return 1;
    // return 2;
    return 3;
    } // P=3
  case 4: // DEN DEBUG: check it carefully and perform calculation of coefficients in long double
    {
    if (N_pairs_of_complex_roots==2) // no real roots
      {
      re=rnr(rng); while(std::abs(im=rnr(rng))<std::abs(re)) {}
      RE=re; IM=im;
      IM=pr_product_difference(RE,RE,-IM,IM); // RE*RE+IM*IM
      RE*=-2.0L; // irreducible quadratic polynomial is (x^2 + re*x + im)
      u=rnr(rng); while(std::abs(v=rnr(rng))<std::abs(u)) {}
      U=u; V=v;
      V=pr_product_difference(U,U,-V,V); // U*U+V*V
      U*=-2.0L; // irreducible quadratic polynomial is (x^2 + u*x + v)
      // multiply both irreducible quadrics
      coefficients[0]=static_cast<fp_t>(IM*V); coefficients[1]=static_cast<fp_t>(pr_product_difference(RE,V,-IM,U)); // RE*V+IM*U;
      coefficients[2]=static_cast<fp_t>(std::fma(RE,U,IM+V)); // IM+RE*U+V
      coefficients[3]=static_cast<fp_t>(RE+U);
      return 0;
      }
    else if (N_pairs_of_complex_roots==1) // two real roots
      {
      re=rnr(rng); while(std::abs(im=rnr(rng))<std::abs(re)) {}
      RE=re; IM=im;
      IM=pr_product_difference(RE,RE,-IM,IM); // RE*RE+IM*IM
      RE*=-2.0L; // irreducible quadratic polynomial is (x^2 + re*x + im); multiply it by the rest
      // 2 real roots follow
      if (N_clustered_roots==2) // 2 clustered roots
        { roots[0]=u=rnr(rng); v=rnc(rng); roots[1]=v=(u>root_mid_sweep ? u-v : u+v); }
      else if (N_multiple_roots==2) // 2 multiple roots
        { roots[1]=roots[0]=u=v=rnr(rng); }
      else // 2 distinct roots
        { roots[0]=u=rnr(rng); roots[1]=v=rnr(rng); }
      U=u; V=v;
      TMP=-U-V; V*=U; U=TMP; // two-real-root quadratic polynomial is (x^2 + u*x + v)
      // multiply irreducible and reducible quadrics
      coefficients[0]=static_cast<fp_t>(IM*V); coefficients[1]=static_cast<fp_t>(pr_product_difference(RE,V,-IM,U)); // RE*V+IM*U
      coefficients[2]=static_cast<fp_t>(std::fma(RE,U,IM+V)); // IM+RE*U+V
      coefficients[3]=static_cast<fp_t>(RE+U);
      return 2;
      }
    else if (N_clustered_roots==4) // 4 clustered roots
      {
      roots[0]=re=rnr(rng); im=rnc(rng); u=rnc(rng); v=rnc(rng);
      roots[1]=im=(re>root_mid_sweep ? (roots[3]=v=(re-im-u-v), roots[2]=u=(re-im-u), re-im) :
                                       (roots[3]=v=(re+im+u+v), roots[2]=u=(re+im+u), re+im) );
      }
    else if (N_clustered_roots==3) // 3 clustered roots and 1 single root
      {
      roots[0]=re=rnr(rng); im=rnc(rng); u=rnc(rng);
      roots[1]=im=(re>root_mid_sweep ? (roots[2]=u=(re-im-u), re-im) : (roots[2]=u=(re+im+u), re+im) );
      roots[3]=v=rnr(rng); // a single root
      }
    else if (N_clustered_roots==2) // 2 clustered roots
      {
      roots[0]=re=rnr(rng); im=rnc(rng); roots[1]=im=(re>root_mid_sweep ? re-im : re+im);
      if (N_multiple_roots==2) // 2 multiple roots
        { roots[3]=roots[2]=v=u=rnr(rng); }
      else // 2 single roots
        { roots[2]=u=rnr(rng); roots[3]=v=rnr(rng); }
      }
    else if (N_multiple_roots==4) // 4 multiple roots
      { roots[3]=roots[2]=roots[1]=roots[0]=v=u=im=re=rnr(rng); }
    else if (N_multiple_roots==3) // 3 multiple roots and 1 single root
      { roots[2]=roots[1]=roots[0]=u=im=re=rnr(rng); roots[3]=v=rnr(rng); }
    else if (N_multiple_roots==2) // 2 multiple roots and 2 single roots
      { roots[1]=roots[0]=im=re=rnr(rng); roots[2]=u=rnr(rng); roots[3]=v=rnr(rng); }
    else // 4 distinct single roots
      { roots[0]=re=rnr(rng); roots[1]=im=rnr(rng); roots[2]=u=rnr(rng); roots[3]=v=rnr(rng); }
    // compute coefficients from 4 roots: re, im, u, v
    RE=re; IM=im; U=u; V=v;
    TMP=-RE-IM; IM*=RE; RE=TMP; // now we have the 1.st quadratic polynomial: x^2 + x*re + im
    TMP=-U-V; V*=U; U=TMP; // now we have the 2.nd quadratic polynomial: x^2 + x*u + v
    coefficients[0]=static_cast<fp_t>(IM*V); coefficients[1]=static_cast<fp_t>(pr_product_difference(RE,V,-IM,U)); // RE*V+IM*U
    coefficients[2]=static_cast<fp_t>(std::fma(RE,U,IM+V)); // IM+RE*U+V
    coefficients[3]=static_cast<fp_t>(RE+U); return 4;
    } // P=4
  } // switch (P)
return -1; // unreachable, means a flaw in control here
}

// Compares two vectors of roots; root orderings play no role. For each entry in (roots_ground_truth),
// the closest entry in (roots_to_check) is found and coresponding distance found. Among such distances
// the largest will be stored to (max_deviation)
template<typename fp_t> int compare_roots(
unsigned N_roots_to_check, // number of roots in (roots_to_check)
unsigned N_roots_ground_truth,  // number of roots in (roots_ground_truth)
std::vector<fp_t> &roots_to_check, // one should take into account only first (N_roots_to_check) roots here
std::vector<fp_t> &roots_ground_truth, // one should take into account only first (N_roots_ground_truth) roots here
fp_t &max_absolute_error, // here the greatest among the smallest deviations of the roots in (roots_to_check) and (roots_ground_truth)
// will be placed
fp_t &max_relative_error){
    int rv = (N_roots_to_check<N_roots_ground_truth) ? PR_AT_LEAST_ONE_ROOT_LOST :
      ( (N_roots_to_check>N_roots_ground_truth) ? PR_AT_LEAST_ONE_ROOT_IS_FAKE : PR_NUMBERS_OF_ROOTS_EQUAL );
    long double abs = std::numeric_limits<long double >::max();
    long double  rel = std::numeric_limits<long double >::max();
    auto size = roots_to_check.size();
    for(int j = 0;j<size; j++)
    for(int i = 0;i < size; i++){
        long double  absLoc = std::abs((long double)(roots_ground_truth[i])-(long double)(roots_to_check[(i + j) % size]));
        abs = std::min(absLoc,abs);
        rel = std::min(std::abs(
                (long double)(absLoc + std::numeric_limits<fp_t>::epsilon())/
                        (long double)(std::max(roots_to_check[(i + j) % size],roots_ground_truth[i]) + std::numeric_limits<fp_t>::epsilon())),rel);
    }
    max_absolute_error = abs;
    max_relative_error = rel;
    return rv;
}

// checks attainable number of real roots in a polynomial: a*x^4 + b*x^3 + c*x^2 + d*x + e; multiple root is treated as separate roots
template <typename fp_t> int number_of_roots(unsigned P, // polynomial degree
fp_t a, fp_t b, fp_t c, fp_t d, fp_t e) // polynomial coefficients
{
switch (P)
  {
  case 1:
    { fp_t x1=-e/d; return std::isinf(x1) ? 0 : 1; break; }
  case 2:
    {
    fp_t x1=d/(static_cast<fp_t>(-2.0L)*c);
    if (std::isinf(x1)) // technically is degenerated quadric
      { x1=-e/d; return std::isinf(x1) ? 0 : 1; }
    fp_t y1=std::fma(c,x1,d); y1=std::fma(y1,x1,e); // y1=c*x1^2+d*x1+e
    return c*y1<static_cast<fp_t>(0.0L) ? 2 : 0;
    break;
    }
  case 3:
    std::cerr << "P=" << " not implemented";
    assert(0);
  }
return -1; // inaccesible
}

typedef float fp_t; // floating point type for all the operations
//typedef double fp_t; // floating point type for all the operations

int main(void)
{
// parameter section //////////////////////////////////////////////////////////////////////////////
unsigned P=3; // max power to test (min is zero)
unsigned long long N_tests=NTESTS; // total number of tests based on makefile build
unsigned N_test_max_to_verbose_output=10; // maximal number of tests for which extended-verbosity output is produced
  // for every polynomial test
unsigned N_pairs_of_complex_roots=0; // how many pairs of complex conjugate roots to introduce in each test
unsigned N_clustered_roots=3; // how many clustered roots to introduce in each test; all the clustered roots are real
unsigned N_multiple_roots=0; // how many multiple roots to introduce in each test; all multiple roots are real
fp_t max_distance_between_clustered_roots=1e-5; // maximal distance between the closest of the clustered roots in each test
  // (efficient if N_clustered_roots>=2)
fp_t root_sweep_low=-1.0, root_sweep_high=1.0; // low and high boundaries of all real roots; imaginary parts of complex conjugate roots
  // are created in the same range
fp_t complex_value_tolerance=0.0; // a very small non-negative real number that sets the tolerance to complex roots having
  // tiny imagonary part and forces them to become multiple real roots; exact root classification is done with complex_value_tolerance=0
  // but this may lead to loosing very close truly real roots due to imprecision of floating point arithmetics
fp_t value_at_root_tolerance=0.0; // a very small non-negative real number that sets the tolerance to the departation from zero of the polynomial value at "double root candidate" point
bool exclude_lost_root_cases_from_statistics=true; // do not count in statistics and analyze the failures to find all the roots
bool test_for_precision=true; // do check correctness of pwr(.)
pr_discriminant_computation_method=PR_DISCRIMINANT_USE_TRADITIONAL_OPERATIONS_NORMALIZED; //PR_DISCRIMINANT_USE_FMA_UNNORMALIZED; // selector of discriminant computation method for quadratic polynomial
pr_lost_roots_saving_1=1; // true allows enabling the 1.st method for the detection of the lost roots (global variable)
pr_lost_roots_saving_2=1; // true allows enabling the 2.nd method for the detection of the lost roots (global variable)

// code section ///////////////////////////////////////////////////////////////////////////////////
std::vector<fp_t> roots_gt_this_test(P), roots_found_this_test(P), coefficients_this_test(P+1), // current polynomial
  roots_gt_all_tests_ae_worst(P), roots_found_all_tests_ae_worst(P), coefficients_all_tests_ae_worst(P+1), // {gt roots, found roots, coefficients} of the worst polynomial for absolute error
  roots_gt_all_tests_re_worst(P), roots_found_all_tests_re_worst(P), coefficients_all_tests_re_worst(P+1); // {gt roots, found roots, coefficients} of the worst polynomial for relative error

std::vector<long double> roots_found_this_test_ld(P);  
fp_t ae_this_test_worst, ae_all_tests_worst=static_cast<fp_t>(-1.0L), // absolute error: an impossible value that will be updated for sure
     re_this_test_worst, re_all_tests_worst=static_cast<fp_t>(-1.0L); // relative error: an impossible value that will be updated for sure
long double poly; // evaluated polynomial value
int i, j, rv, N_roots_gt_this_test, N_roots_found_this_test, N_roots_found_this_test_ld, N_roots_verified_this_test, N_roots_verified_this_test_ld,
  N_roots_gt_all_tests_ae_worst=0, N_roots_found_all_tests_ae_worst=0, // {ground truth, found} number of roots of the worst polynomial for absolute error
  N_roots_gt_all_tests_re_worst=0, N_roots_found_all_tests_re_worst=0, // {ground truth, found} number of roots of the worst polynomial for relative error
  N_true_roots_lost=0, N_fake_roots_added=0; // total counters of {lost, fake} roots aover all tests
unsigned long long n_test;
N_additive_feedback_fired=0ULL; N_multiplicative_feedback_fired=0ULL;

// Baydoun or Vieta, based on make build
METHOD;

if (test_for_precision) // check correctness
  {
  std::streamsize fp_precision_original=std::cout.precision(); // save default precision to provide maximal reasonable output
  for (n_test=0; n_test<N_tests; ++n_test) // N_tests
    {
    roots_found_this_test.clear();
    std::vector<std::complex<fp_t>> b_roots;
    N_roots_gt_this_test=generate_polynomial<fp_t>(P, N_pairs_of_complex_roots, N_clustered_roots, N_multiple_roots,
      max_distance_between_clustered_roots, root_sweep_low, root_sweep_high, roots_gt_this_test, coefficients_this_test);

    
    Solver(coefficients_this_test, b_roots, true);
    std::for_each(b_roots.begin(), b_roots.end(), [&roots_found_this_test](std::complex<fp_t> x){
        // if(x.imag() <= 1){
          // std::cout << "dddd " << x.real();
          roots_found_this_test.push_back(x.real());
        // }
    });
    N_roots_found_this_test = roots_found_this_test.size();

    rv=compare_roots<fp_t>(N_roots_found_this_test, N_roots_gt_this_test, roots_found_this_test, roots_gt_this_test,
      ae_this_test_worst, re_this_test_worst);

    N_roots_verified_this_test=P; /* number_of_roots<fp_t>(P,coefficients_this_test[4],coefficients_this_test[3],
      coefficients_this_test[2],coefficients_this_test[1],coefficients_this_test[0]); DEBUG */
    N_roots_verified_this_test_ld=P; /* number_of_roots<long double>(P,static_cast<long double>(coefficients_this_test[4]),
      static_cast<long double>(coefficients_this_test[3]),static_cast<long double>(coefficients_this_test[2]),
      static_cast<long double>(coefficients_this_test[1]),static_cast<long double>(coefficients_this_test[0])); DEBUG */

    if (N_tests<=N_test_max_to_verbose_output)
      {
      std::cout << std::endl << "P=" << P << ", test No " << n_test << " out of " << N_tests << std::endl
        << "comparison return value=" << (rv==PR_NUMBERS_OF_ROOTS_EQUAL ? "PR_NUMBERS_OF_ROOTS_EQUAL" :
          (rv==PR_AT_LEAST_ONE_ROOT_LOST ? "PR_AT_LEAST_ONE_ROOT_LOST" : 
          (rv==PR_AT_LEAST_ONE_ROOT_IS_FAKE ? "PR_AT_LEAST_ONE_ROOT_IS_FAKE" : "UNKNOWN VALUE")))
        << ", N_roots_gt_this_test=" << N_roots_gt_this_test << ", N_roots_found_this_test=" << N_roots_found_this_test
        << ", N_roots_verified_this_test=" << N_roots_verified_this_test
        << ", N_roots_verified_this_test_ld=" << N_roots_verified_this_test_ld << std::endl
        << "ae_this_test_worst=" << ae_this_test_worst
        << ", re_this_test_worst=" << re_this_test_worst << " ------------------------------" << std::endl;

      for (i=0; i<N_roots_gt_this_test; ++i)
        {
        poly=std::fma<long double>(coefficients_this_test[P],roots_gt_this_test[i],coefficients_this_test[P-1]);
        for (j=P-2; j>=0; --j) poly=std::fma<long double>(poly,roots_gt_this_test[i],coefficients_this_test[j]);
        std::cout << std::setprecision(std::numeric_limits<fp_t>::digits10 + 3)
          << "p(r_gt[" << i << "]=" << roots_gt_this_test[i] << ")=" << poly << (i<N_roots_gt_this_test-1 ? ", " : "\n");
        }
      for (i=0; i<N_roots_found_this_test; ++i)
        {
        poly=std::fma<long double>(coefficients_this_test[P],roots_found_this_test[i],coefficients_this_test[P-1]);
        for (j=P-2; j>=0; --j) poly=std::fma<long double>(poly,roots_found_this_test[i],coefficients_this_test[j]);
        std::cout << std::setprecision(std::numeric_limits<fp_t>::digits10 + 3)
          << "p(r_f[" << i << "]=" << roots_found_this_test[i] << ")=" << poly << (i<=N_roots_found_this_test-1 ? ", " : "\n");
        }
      for (i=0; i<=P; ++i) std::cout << std::setprecision(std::numeric_limits<fp_t>::digits10 + 3)
        << "c[" << i << "]=" << coefficients_this_test[i] << (i<P ? ", " : "\n");
      }

    N_true_roots_lost+=(N_roots_found_this_test<N_roots_gt_this_test)*(N_roots_gt_this_test-N_roots_found_this_test);
    N_fake_roots_added+=(N_roots_found_this_test>N_roots_gt_this_test)*(N_roots_found_this_test-N_roots_gt_this_test);
    if (ae_this_test_worst>ae_all_tests_worst && (!exclude_lost_root_cases_from_statistics || N_roots_found_this_test>=N_roots_gt_this_test))
      {
      roots_gt_all_tests_ae_worst=roots_gt_this_test; N_roots_gt_all_tests_ae_worst=N_roots_gt_this_test;
      roots_found_all_tests_ae_worst=roots_found_this_test; N_roots_found_all_tests_ae_worst=N_roots_found_this_test;
      coefficients_all_tests_ae_worst=coefficients_this_test; ae_all_tests_worst=ae_this_test_worst;
      }
    if (re_this_test_worst>re_all_tests_worst && (!exclude_lost_root_cases_from_statistics || N_roots_found_this_test>=N_roots_gt_this_test))
      {
      roots_gt_all_tests_re_worst=roots_gt_this_test; N_roots_gt_all_tests_re_worst=N_roots_gt_this_test;
      roots_found_all_tests_re_worst=roots_found_this_test; N_roots_found_all_tests_re_worst=N_roots_found_this_test;
      coefficients_all_tests_re_worst=coefficients_this_test; re_all_tests_worst=re_this_test_worst;
      }

    } // for (n_test=0; n_test<N_tests; ++n_test)

  // output main settings
  std::cout << std::endl << "--------------- settings ----------------" << std::endl
    << "P=" << P << ", N_tests=" << N_tests << ", N_pairs_of_complex_roots=" << N_pairs_of_complex_roots
    << ", N_clustered_roots=" << N_clustered_roots << ", N_multiple_roots=" << N_multiple_roots << std::endl
    << "max_distance_between_clustered_roots=" << max_distance_between_clustered_roots
    << ", root_sweep=[" << root_sweep_low << "..." << root_sweep_high << "]\ncomplex_value_tolerance=" << complex_value_tolerance
    << ", value_at_root_tolerance=" << value_at_root_tolerance << ", exclude_lost_root_cases_from_statistics=" << exclude_lost_root_cases_from_statistics << std::endl
    << "pr_lost_roots_saving_1=" << pr_lost_roots_saving_1 << ", pr_lost_roots_saving_2=" << pr_lost_roots_saving_2;
  if (P==2)
    {
    std::cout << ", pr_discriminant_computation_method="
      << (::pr_discriminant_computation_method==PR_DISCRIMINANT_USE_TRADITIONAL_OPERATIONS_NORMALIZED ? "PR_DISCRIMINANT_USE_TRADITIONAL_OPERATIONS_NORMALIZED" :
         (::pr_discriminant_computation_method==PR_DISCRIMINANT_USE_FMA_UNNORMALIZED ? "PR_DISCRIMINANT_USE_FMA_UNNORMALIZED" :
         (::pr_discriminant_computation_method==PR_DISCRIMINANT_USE_FMA_NORMALIZED ? "PR_DISCRIMINANT_USE_FMA_NORMALIZED" : "UNKNOWN")));
    }

  // output total statistics
  std::cout << std::endl << std::endl << "--------------- statistics ----------------" << std::endl
    << "N_true_roots_lost=" << N_true_roots_lost << ", N_fake_roots_added=" << N_fake_roots_added
      << ", N_additive_feedback_fired=" << N_additive_feedback_fired
      << ", N_multiplicative_feedback_fired=" << N_multiplicative_feedback_fired << std::endl;

  // output the worst cases
  std::cout << std::endl << "--------------- worst case: absolute error: ----------------" << std::endl
    << "ae_all_tests_worst=" << ae_all_tests_worst
    << ", N_roots_gt_all_tests_ae_worst=" << N_roots_gt_all_tests_ae_worst
    << ", N_roots_found_all_tests_ae_worst=" << N_roots_found_all_tests_ae_worst << std::endl;
  /*
  for (i=0; i<N_roots_gt_all_tests_ae_worst; ++i) std::cout << "r_gt_ae[" << i << "]="
    << std::setprecision(std::numeric_limits<fp_t>::digits10 + 3) << roots_gt_all_tests_ae_worst[i]
    << (i<N_roots_gt_all_tests_ae_worst-1 ? ", " : "\n");
  for (i=0; i<N_roots_found_all_tests_ae_worst; ++i) std::cout << "r_f_ae[" << i << "]="
    << std::setprecision(std::numeric_limits<fp_t>::digits10 + 3) << roots_found_all_tests_ae_worst[i]
    << (i<N_roots_found_all_tests_ae_worst-1 ? ", " : "\n");
  for (i=0; i<=P; ++i) std::cout << std::setprecision(std::numeric_limits<fp_t>::digits10 + 1)
    << "c_ae[" << i << "]=" << coefficients_all_tests_ae_worst[i] << (i<P ? ", " : "\n");
  */

  std::cout << std::endl << "--------------- worst case: relative error: ----------------" << std::endl
    << "re_all_tests_worst=" << re_all_tests_worst
    << ", N_roots_gt_all_tests_re_worst=" << N_roots_gt_all_tests_re_worst
    << ", N_roots_found_all_tests_re_worst=" << N_roots_found_all_tests_re_worst << std::endl;
  /*
  for (i=0; i<N_roots_gt_all_tests_re_worst; ++i) std::cout << "r_gt_re[" << i << "]="
    << std::setprecision(std::numeric_limits<fp_t>::digits10 + 3) << roots_gt_all_tests_re_worst[i]
    << (i<N_roots_gt_all_tests_re_worst ? ", " : "\n");
  for (i=0; i<N_roots_found_all_tests_re_worst; ++i) std::cout << "r_f_re[" << i << "]="
    << std::setprecision(std::numeric_limits<fp_t>::digits10 + 3) << roots_found_all_tests_re_worst[i]
    << (i<N_roots_found_all_tests_re_worst ? ", " : "\n");
  for (i=0; i<=P; ++i) std::cout << std::setprecision(std::numeric_limits<fp_t>::digits10 + 1)
    << "c_re[" << i << "]=" << coefficients_all_tests_re_worst[i] << (i<P ? ", " : "\n");
  */

  std::cout << "\nWorst case method results ae: ";
  for (auto &root : roots_found_all_tests_ae_worst){
      std::cout << root << " ";
  }
  std::cout << "\nGround truth worst ae: {";
  for (auto &root : roots_gt_all_tests_ae_worst){
      std::cout << root << ", ";
  }
  std::cout << "}\nCoeffs worst ae: {";
  for (auto &root : coefficients_all_tests_ae_worst){
      std::cout << root << ", ";
  }


  std::cout << "}\nWorst case method results re: ";
  for (auto &root : roots_found_all_tests_re_worst){
      std::cout << root << " ";
  }
  std::cout << "\nGround truth worst re: {";
  for (auto &root : roots_gt_all_tests_re_worst){
      std::cout << root << ", ";
  }
  std::cout << "}\nCoeffs worst re: {";
  for (auto &root : coefficients_all_tests_re_worst){
      std::cout << root << ", ";
  }
  std::cout << "}\n";

  std::cout << std::setprecision(fp_precision_original); // restore default precision
  } // if (test_for_precision)

return 0;
}
#endif
