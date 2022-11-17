#include <iostream>
#include "excerpt.h"
#include "Methods.h"

int main(int argc,char* argv[]){
    std::setprecision(30);
    int count;
    if (argc == 1){
        std::cout << "Введите количество полиномов > ";
        std::cin >> count;
    }
    else{
        count = std::atoi(argv[1]);
    }
    std::cout << "Генерирую полиномы.\n";

    Baydoun<double> SolverB;
    Vieta<double> SolverV;
    for(int i = 0; i < count; i++){
        std::vector<double> roots(3), coefficients(4);
        generate_polynomial<double>(3, 0, 3, 0, 1e-5, -1, 1, roots, coefficients);
        std::cout << "Коэффициенты:\n";
        std::cout << "x^3 + (" << coefficients[2] << ")x^2 + (" << coefficients[1] <<
            ")x + (" << coefficients[0] << ") = 0\n";
        std::cout << "Корни\n";
        for (auto &root : roots){
            std::cout << "(x" << ((root > 0) ? "-" : "+") << fabs(root) << ")";
        }
        std::cout << "Результаты\n";
        std::vector<std::complex<double>> b_roots;
        int roots_num = SolverB(coefficients, b_roots, true);
        std::cout << "Baydoun. Найдено " << roots_num << " корня.\n";
        for (auto &root : b_roots){
            std::cout << root << "\n";
        }
        std::vector<std::complex<double>> v_roots;
        roots_num = SolverV(coefficients, v_roots, true);
        std::cout << "Vieta. Найдено " << roots_num << " корня.\n";
        for (auto &root : v_roots){
            std::cout << root << "\n";
        }
        std::cout << "\n==================================================\n";
    }
    
}
