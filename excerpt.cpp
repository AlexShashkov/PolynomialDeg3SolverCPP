#include "excerpt.h"

int pr_discriminant_computation_method = 0;
int pr_lost_roots_saving_1 = 0;
int pr_lost_roots_saving_2 = 0;
unsigned long long N_additive_feedback_fired = 0ULL;
unsigned long long N_multiplicative_feedback_fired = 0ULL;


/* computes (a*b - c*d) with precision not worse than 1.5*(unit of the least precision) suggested in Claude-Pierre Jeannerod,
Nicolas Louvet, and Jean-Michel Muller, "Further Analysis of Kahan's Algorithm for the Accurate Computation of 2x2 Determinants".
Mathematics of Computation, Vol. 82, No. 284, Oct. 2013, pp. 2245-2264 */
template<typename fp_t>
inline fp_t pr_product_difference(fp_t a, fp_t b, fp_t c, fp_t d) {
    auto tmp = d * c;
    return std::fma(a, b, -tmp) + std::fma(-d, c, tmp);
}


// Finds all real roots of a quadratic polynomial in the form: a*x^2 + b*x + c = 0
// This implementation should process correctly all the situations and returns 0,1 (linear case, a=0 !) or 2 roots. It returns
// double quadratic root twice, in order to keep statistics more correct way; it serves solely verification purposes because
// no optimizations and enhancements introduced
int polyroots2plain(float a, float b, float c, // polynomial coefficients
                    std::vector<float> &roots) // pre allocated vector where to put the roots; the size should not be less than 2
{
    auto tmp = 0.0F, D = 0.0F;
    if (a != 0.0F) {
        tmp = b / a;
        D = c / a;
    }
    if (a == 0.0F || std::isinf(std::abs(tmp)) ||
        std::isinf(std::abs(D))) // exactly or technically, the polynomial is not quadratic: believe a=0
    {
        if (b == 0.0F) return 0;
        tmp = -c / b;
        if (std::isinf(std::abs(tmp))) return 0;
        roots[0] = tmp;
        return 1;
    }

// non-degenerate quadratic polynomial is resulted and further on processed
    D = b * b - 4.0F * a * c; // discriminant

    if (D >= 0.0F) // 2 distinct real roots or a single root
    {
        roots[0] = tmp = (-b - std::copysignf(std::sqrt(D), b)) / (2.0F * a);
        if (D == 0.0F) {
            roots[1] = tmp;
            return 2;
        }
        roots[1] = c / tmp;
        if (std::isinf(std::abs(roots[1]))) roots[1] = (-b + std::copysignf(std::sqrt(D), b)) / (2.0F * a);
        return 2;
    } else return 0; // D<0.0
}


// Creates a test polynomial, both in the form of roots, e.g. (x-roots[0])*(x-roots[1])*(quadratic polynomial with no real roots)
// and represented by its coefficients, e.g.
// (coefficients[4]=1)*x^4 + coefficients[3]*x^3 + coefficients[2]*x^2 + coefficients[1]*x + coefficients[0].
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
        std::vector<fp_t> &coefficients) // storage where to put the coefficients; size should exceed P
{
    int n_simple_roots = P - 2 * N_pairs_of_complex_roots - N_clustered_roots - N_multiple_roots;
    assert(N_clustered_roots != 1);
    assert(N_multiple_roots != 1);
    assert(n_simple_roots >= 0 && P <= 4);
    assert(root_sweep_high - root_sweep_low > 2 * P * max_distance_between_clustered_roots);

    unsigned long long seed =
            std::chrono::system_clock::now().time_since_epoch().count() + std::rand(); // counts milliseconds
    std::mt19937_64 rng(seed); // randomize seed from the clock
    std::uniform_real_distribution<fp_t> rnr(root_sweep_low,
                                             root_sweep_high); // uniform random data generator for single roots
    std::uniform_real_distribution<fp_t> rnc(static_cast<fp_t>(0.0L),
                                             max_distance_between_clustered_roots); // uniform random data generator for root clusters
    fp_t re, im, u, v, root_mid_sweep = root_sweep_low + 0.5 * (root_sweep_high - root_sweep_low);
    long double RE, IM, U, V; // high-precisioon counterparts of re, im, u, v

    coefficients[P] = static_cast<fp_t>(1.0L); // invariant
    switch (P) {
        case 0:
            coefficients[0] = rnr(rng);
            return 0;
        case 1:
            coefficients[0] = -(roots[0] = rnr(rng));
            return 1;
        case 2: {
            if (N_pairs_of_complex_roots == 1) // no real roots
            {
                re = rnr(rng);
                while ((im = rnr(rng)) == static_cast<fp_t>(0.0L)) {}
                RE = re;
                IM = im;
                coefficients[1] = static_cast<fp_t>(-2.0L * RE); // -2*re
                coefficients[0] = static_cast<fp_t>(pr_product_difference(RE, RE, -IM, IM)); // re*re+im*im
                return 0;
            } else if (N_clustered_roots == 2) // 2 close but distinct roots
            {
                roots[0] = re = rnr(rng);
                while ((im = rnc(rng)) == static_cast<fp_t>(0.0L)) {}
                roots[1] = im = (re >= root_mid_sweep ? re - im : re + im);
            } else if (N_multiple_roots == 2) // double root counted as a single root
            { roots[1] = roots[0] = im = re = rnr(rng); }
            else // 2 distinct single roots
            {
                roots[0] = re = rnr(rng);
                while ((im = rnr(rng)) == re) {}
                roots[1] = im = rnr(rng);
            }
            RE = re;
            IM = im;

            coefficients[1] = static_cast<fp_t>( -RE - IM );
            coefficients[0] = static_cast<fp_t>( RE * IM );
            return 2; // return ((re!=im) ? 2 : 1);
        } // P=2
        case 3: {
            if (N_pairs_of_complex_roots == 1) // one real root
            {
                re = rnr(rng);
                while ((im = rnr(rng)) == static_cast<fp_t>(0.0L)) {}
                roots[0] = u = rnr(rng);
                RE = re;
                IM = im;
                U = u;

                IM = pr_product_difference(RE, RE, -IM, IM); // re*re+im*im
                RE *= -2.0L; // irreducible quadratic polynomial is (x^2 + re*x + im); multiply it by (x-u)
                coefficients[0] = static_cast<fp_t>(-IM * U);
                coefficients[2] = static_cast<fp_t>(RE - U);
                coefficients[1] = static_cast<fp_t>(std::fma(-RE, U, IM)); // im-re*u;
                return 1;
            } else if (N_clustered_roots == 3) // 3 clustered distinct roots
            {
                roots[0] = re = rnr(rng);
                while ((im = rnc(rng)) == static_cast<fp_t>(0.0L)) {}
                while ((u = rnc(rng)) == static_cast<fp_t>(0.0L)) {}
                roots[1] = im = (re > root_mid_sweep ? roots[2] = u = (re - im - u), re - im : roots[2] = u = (re + im +
                                                                                                               u), re +
                                                                                                                   im);
            } else if (N_clustered_roots == 2) // 2 clustered roots, 1 single root; all distinct
            {
                roots[0] = re = rnr(rng);
                while ((im = rnc(rng)) == static_cast<fp_t>(0.0L)) {}
                roots[1] = im = (re > root_mid_sweep ? re - im : re + im);
                do { roots[2] = u = rnr(rng); } while (u == re || u == roots[1]);
            } else if (N_multiple_roots == 3) // triple root counted as a single root
            { roots[2] = roots[1] = roots[0] = u = im = re = rnr(rng); }
            else if (N_multiple_roots == 2) // double root and 1 single root; totally 2 roots
            {
                roots[1] = roots[0] = im = re = rnr(rng);
                while ((roots[2] = u = rnr(rng)) == re) {}
            } else // 3 distinct single roots
            {
                roots[0] = re = rnr(rng);
                while ((roots[1] = im = rnr(rng)) == re) {}
                do { roots[2] = u = rnr(rng); } while (u == re || u == im);
            }
            RE = re;
            IM = im;
            U = u;
            coefficients[2] = static_cast<fp_t>(-RE - IM - U);
            coefficients[0] = static_cast<fp_t>(-RE * IM * U);
            V = pr_product_difference(RE, IM, -RE, U);
            coefficients[1] = static_cast<fp_t>(std::fma(IM, U, V)); // re*im+re*u+im*u=im*u+(re*im-(-re*u));
            // if (re!=im && im!=u && u!=re) return 3;
            // if (re==im && im==u) return 1;
            // return 2;
            return 3;
        } // P=3
        case 4: // DEN DEBUG: check it carefully and perform calculation of coefficients in long double
        {
            if (N_pairs_of_complex_roots == 2) // no real roots
            {
                re = rnr(rng);
                im = rnr(rng);
                im = re * re + im * im;
                re *= static_cast<fp_t>(-2.0L); // irreducible quadratic polynomial is (x^2 + re*x + im)
                u = rnr(rng);
                v = rnr(rng);
                v = u * u + v * v;
                u *= static_cast<fp_t>(-2.0L); // irreducible quadratic polynomial is (x^2 + u*x + v)
                // multiply both irreducible quadratics
                coefficients[0] = im * v;
                coefficients[1] = re * v + im * u;
                coefficients[2] = im + re * u + v;
                coefficients[3] = re + u;
                return 0;
            } else if (N_pairs_of_complex_roots == 1) // two real roots
            {
                re = rnr(rng);
                im = rnr(rng);
                roots[0] = u = rnr(rng);
                im = re * re + im * im;
                re *= static_cast<fp_t>(-2.0L); // irreducible quadratic polynomial is (x^2 + re*x + im); multiply it by (x-u)

                if (N_clustered_roots == 2) // 2 clustered roots
                {
                    roots[0] = u = rnr(rng);
                    v = rnc(rng);
                    roots[1] = v = (u > root_mid_sweep ? u - v : u + v);
                } else if (N_multiple_roots == 2) // 2 multiple roots
                { roots[1] = roots[0] = u = v = rnr(rng); }
                else // 2 distinct roots
                {
                    roots[0] = u = rnr(rng);
                    roots[1] = v = rnr(rng);
                }
                root_mid_sweep = -u - v;
                v *= u;
                u = root_mid_sweep; // two-root quadratic polynomial is (x^2 + u*x + v)
                // multiply irreducible and reducible quadrics
                coefficients[0] = im * v;
                coefficients[1] = re * v + im * u;
                coefficients[2] = im + re * u + v;
                coefficients[3] = re + u;
                return 2;
            } else if (N_clustered_roots == 4) // 4 clustered roots
            {
                roots[0] = re = rnr(rng);
                im = rnc(rng);
                u = rnc(rng);
                v = rnc(rng);
                roots[1] = im = (re > root_mid_sweep ? roots[3] = v = (re - im - u - v), roots[2] = u = (re - im - u),
                        re - im :
                        roots[3] = v = (re + im + u + v), roots[2] = u = (re + im + u), re + im);
            } else if (N_clustered_roots == 3) // 3 clustered roots and 1 single root
            {
                roots[0] = re = rnr(rng);
                im = rnc(rng);
                u = rnc(rng);
                roots[1] = im = (re > root_mid_sweep ? roots[2] = u = (re - im - u), re - im : roots[2] = u = (re + im +
                                                                                                               u), re +
                                                                                                                   im);
                roots[3] = v = rnc(rng); // 1 single root
            } else if (N_clustered_roots == 2) // 2 clustered roots
            {
                roots[0] = re = rnr(rng);
                im = rnc(rng);
                roots[1] = im = (re > root_mid_sweep ? re - im : re + im);
                if (N_multiple_roots == 2) // 2 multiple roots
                { roots[3] = roots[2] = v = u = rnr(rng); }
                else // 2 single roots
                {
                    roots[2] = u = rnr(rng);
                    roots[3] = v = rnr(rng);
                }
            } else if (N_multiple_roots == 4) // 4 multiple roots
            { roots[3] = roots[2] = roots[1] = roots[0] = v = u = im = re = rnr(rng); }
            else if (N_multiple_roots == 3) // 3 multiple roots and 1 single root
            {
                roots[2] = roots[1] = roots[0] = u = im = re = rnr(rng);
                roots[3] = v = rnr(rng);
            } else if (N_multiple_roots == 2) // 2 multiple roots and 2 single roots
            {
                roots[1] = roots[0] = im = re = rnr(rng);
                roots[2] = u = rnr(rng);
                roots[3] = v = rnr(rng);
            } else // 4 distinct single roots
            {
                roots[0] = re = rnr(rng);
                roots[1] = im = rnr(rng);
                roots[2] = u = rnr(rng);
                roots[3] = v = rnr(rng);
            }
            // compute coefficients from 4 roots: re, im, u, v
            root_mid_sweep = -re - im;
            im *= re;
            re = root_mid_sweep; // now we have the 1.st quadratic polynomial: x^2 + x*re + im
            root_mid_sweep = -u - v;
            v *= u;
            u = root_mid_sweep; // now we have the 2.nd quadratic polynomial: x^2 + x*u + v
            coefficients[0] = im * v;
            coefficients[1] = re * v + im * u;
            coefficients[2] = re * u + im + v;
            coefficients[3] = re + u;
            return 4;
        } // P=4
        default:
            break;
    } // switch (P)
    return -1; // unreachable, means a flaw in control here
}


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
        fp_t &max_relative_error) // here the greatest relative error among all the roots found will be placed
{
    auto absolute_error_max = static_cast<fp_t>(0.0L), relative_error_max = static_cast<fp_t>(0.0L);
    auto rv = (N_roots_to_check < N_roots_ground_truth) ? PR_AT_LEAST_ONE_ROOT_LOST :
              ((N_roots_to_check > N_roots_ground_truth) ? PR_AT_LEAST_ONE_ROOT_IS_FAKE : PR_NUMBERS_OF_ROOTS_EQUAL);
// find the largest distance between the closest pairs of roots: one - from ground truth, one - from found ones
    for (auto i = 0; i < N_roots_ground_truth; ++i) {
        // find the closest found root to the given ground truth root
        auto deviation_min_for_this_root = std::numeric_limits<fp_t>::infinity();
        auto i_closest_root = -1, j_closest_root = -1;
        for (auto j = 0; j < N_roots_to_check; ++j) {
            auto deviation = std::abs(roots_ground_truth[i] - roots_to_check[j]);
            deviation_min_for_this_root =
                    deviation < deviation_min_for_this_root ? i_closest_root = i, j_closest_root = j, deviation
                                                            : deviation_min_for_this_root;
        }
        assert(i_closest_root != -1 and j_closest_root != -1);
        auto relative_error_for_this_root = static_cast<fp_t>(2.0L) * deviation_min_for_this_root /
                                            (std::abs(roots_ground_truth[i_closest_root]) +
                                             std::abs(roots_to_check[j_closest_root]));

        absolute_error_max =
                deviation_min_for_this_root > absolute_error_max ? deviation_min_for_this_root : absolute_error_max;
        relative_error_max =
                relative_error_for_this_root > relative_error_max ? relative_error_for_this_root : relative_error_max;
    }
    max_absolute_error = absolute_error_max;
    max_relative_error = relative_error_max;
    return rv;
}


// checks attainable number of real roots in a polynomial: a*x^4 + b*x^3 + c*x^2 + d*x + e; multiple root is treated as separate roots
template<typename fp_t>
int number_of_roots(unsigned P, // polynomial degree
                    fp_t a, fp_t b, fp_t c, fp_t d, fp_t e) // polynomial coefficients
{
    switch (P) {
        case 1: {
            auto x1 = -e / d;
            return std::isinf(x1) ? 0 : 1;
        }
        case 2: {
            auto x1 = d / (static_cast<fp_t>(-2.0L) * c);
            if (std::isinf(x1)) // technically is degenerated quadratic
            {
                x1 = -e / d;
                return std::isinf(x1) ? 0 : 1;
            }
            auto y1 = std::fma(c, x1, d);
            y1 = std::fma(y1, x1, e); // y1=c*x1^2+d*x1+e
            return c * y1 < static_cast<fp_t>(0.0L) ? 2 : 0;
        }
        case 3:
            std::cerr << "P=" << " not implemented";
            assert(0);
        default:
            break;
    }
    return -1; // inaccessible
}

template<typename fp_t>
int compare_roots_complex(unsigned N_roots_to_check, // number of roots in roots_to_check
                          unsigned N_roots_ground_truth,  // number of roots in roots_ground_truth
                          std::vector<std::complex<fp_t>> &roots_to_check, // one should take into account only first (N_roots_to_check) rots
                          std::vector<fp_t> &roots_ground_truth, // one should take into account only first (N_roots_ground_truth) rots
                          fp_t &max_deviation) {
    std::vector<fp_t> roots_to_check_parsed;
    for (auto root: roots_to_check) {
        if (std::numeric_limits<fp_t>::epsilon() > abs(root.imag())) {
            roots_to_check_parsed.push_back(root.real());
        }
    }
    return compare_roots(roots_to_check_parsed.size(), N_roots_ground_truth, roots_to_check_parsed, roots_ground_truth,
                         max_deviation);
}

template int generate_polynomial<float>(unsigned P, unsigned N_pairs_of_complex_roots, unsigned N_clustered_roots,
                                        unsigned N_multiple_roots, float max_distance_between_clustered_roots,
                                        float root_sweep_low, float root_sweep_high, std::vector<float> &roots,
                                        std::vector<float> &coefficients);

template int generate_polynomial<double>(unsigned P, unsigned N_pairs_of_complex_roots, unsigned N_clustered_roots,
                                         unsigned N_multiple_roots, double max_distance_between_clustered_roots,
                                         double root_sweep_low, double root_sweep_high, std::vector<double> &roots,
                                         std::vector<double> &coefficients);

template int generate_polynomial<long double>(unsigned P, unsigned N_pairs_of_complex_roots, unsigned N_clustered_roots,
                                              unsigned N_multiple_roots,
                                              long double max_distance_between_clustered_roots,
                                              long double root_sweep_low, long double root_sweep_high,
                                              std::vector<long double> &roots,
                                              std::vector<long double> &coefficients);

template int compare_roots<float>(
        unsigned N_roots_to_check,
        unsigned N_roots_ground_truth,
        std::vector<float> &roots_to_check,
        std::vector<float> &roots_ground_truth,
        float &max_absolute_error,
        float &max_relative_error);

template int compare_roots<double>(
        unsigned N_roots_to_check,
        unsigned N_roots_ground_truth,
        std::vector<double> &roots_to_check,
        std::vector<double> &roots_ground_truth,
        double &max_absolute_error,
        double &max_relative_error);

template int compare_roots<long double>(
        unsigned N_roots_to_check,
        unsigned N_roots_ground_truth,
        std::vector<long double> &roots_to_check,
        std::vector<long double> &roots_ground_truth,
        long double &max_absolute_error,
        long double &max_relative_error);
