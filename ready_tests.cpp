#include <iostream>
#include <tuple>

#include "excerpt.h"
#include "Methods.h"

#define number float
using namespace implementations;

std::tuple<number, number> test_method(std::vector<number> &roots, std::vector<number> &real_roots, std::string methodname){
    number max_absolute_error, max_relative_error;
    auto result = compare_roots<number>(roots.size(), real_roots.size(), roots, real_roots, max_absolute_error,
                                max_relative_error);
    switch (result) {
        case PR_2_INFINITE_ROOTS:
            std::cout << "Solution of " << methodname << " contains infinite roots. ";
            break;
        case PR_AT_LEAST_ONE_ROOT_IS_FAKE:
            std::cout << "Solution of " << methodname << " contains at least one fake root. ";
            break;
        case PR_AT_LEAST_ONE_ROOT_LOST:
            std::cout << "At least one root was lost in " << methodname << ". ";
            break;
        default:
            break;
    }

    return {max_absolute_error, max_relative_error};
}

int main(int argc,char* argv[]){
    std::setprecision(30);
    int count;
    number count_as_real = 1;
    number allowed_deviation = 0;

    switch(argc){
        double _car;
        double _ad;
        case 1:
            std::cout << "Input count of polynomials > ";
            std::cin >> count;
            break;
        case 2:
            count = std::atoi(argv[1]);
            break;
        case 3:
            count = std::atoi(argv[1]);
            _car = std::atof(argv[2]);
            count_as_real = static_cast<number>(_car);
            break;
        case 4:
            count = std::atoi(argv[1]);
            _car = std::atof(argv[2]);
            _ad = std::atof(argv[3]);
            count_as_real = static_cast<number>(_car);
            allowed_deviation = static_cast<number>(_ad);
            break;
        default: return 1;
    }

    Baydoun<number> SolverB;
    Vieta<number> SolverV;
    std::vector<number> deviations_abs_B;
    std::vector<number> deviations_abs_V;
    std::vector<number> deviations_rel_B;
    std::vector<number> deviations_rel_V;

    for(int i = 0; i < count; i++){
        std::vector<number> roots(3), coefficients(4);
        generate_polynomial<number>(3, 0, 3, 0, 1e-5, -1, 1, roots, coefficients);
        std::cout << "Coefficients:\n";
        std::cout << "x^3 + (" << coefficients[2] << ")x^2 + (" << coefficients[1] <<
            ")x + (" << coefficients[0] << ") = 0\n";
        std::cout << "True roots\n";
        for (auto &root : roots){
            std::cout << "(x" << ((root > 0) ? "-" : "+") << fabs(root) << ")";
        }
        std::cout << "\nResults\n";
        std::vector<std::complex<number>> b_roots;
        std::vector<number> b_roots_real;

        SolverB(coefficients, b_roots, true);
        std::for_each(b_roots.begin(), b_roots.end(), [count_as_real, &b_roots_real](std::complex<number> x){
            if(x.imag() <= count_as_real) b_roots_real.push_back(x.real());
        });
        std::cout << "Baydoun. Found " << b_roots_real.size() << " roots.\n";
        for (auto &root : b_roots_real){
            std::cout << root << "\n";
        }

        std::vector<std::complex<number>> v_roots;
        std::vector<number> v_roots_real;
        SolverV(coefficients, v_roots, true);
        std::for_each(v_roots.begin(), v_roots.end(), [count_as_real, &v_roots_real](std::complex<number> x){
            if(x.imag() <= count_as_real) v_roots_real.push_back(x.real());
        });
        std::cout << "Vieta. Found " << v_roots_real.size() << " roots.\n";
        for (auto &root : v_roots_real){
            std::cout << root << "\n";
        }
        
        auto [max_abs_B, max_rel_B] = test_method(b_roots_real, roots, "Baydoun");
        auto [max_abs_V, max_rel_V] = test_method(v_roots_real, roots, "Vieta");
        if(max_abs_B > allowed_deviation) deviations_abs_B.push_back(max_abs_B);
        if(max_abs_V > allowed_deviation) deviations_abs_V.push_back(max_abs_V);
        if(max_rel_B > allowed_deviation) deviations_rel_B.push_back(max_rel_B);
        if(max_rel_V > allowed_deviation) deviations_rel_V.push_back(max_rel_V);
        std::cout << "\n==================================================\n";
    }

    auto _max_abs_B = *std::max_element(deviations_abs_B.begin(), deviations_abs_B.end());
    auto _max_abs_V = *std::max_element(deviations_abs_V.begin(), deviations_abs_V.end());
    auto _max_rel_B = *std::max_element(deviations_rel_B.begin(), deviations_rel_B.end());
    auto _max_rel_V = *std::max_element(deviations_rel_V.begin(), deviations_rel_V.end());

    std::cout << "Biggest absolute deviation in Baydoun: " << _max_abs_B << "\n";
    std::cout << "Biggest relative deviation in Baydoun: " << _max_abs_V << "\n";
    std::cout << "Biggest absolute deviation in Vieta: "   << _max_rel_B << "\n";
    std::cout << "Biggest relative deviation in Vieta: "   << _max_rel_V << "\n";
}
