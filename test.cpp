#include "Methods.h"
#include <vector>
#include <iomanip>

#include <random>

/*
double* generateEpsilon(double from, double to, double epsilon, double shakefrom, double shaketo, std::mt19937 gen, std::vector<std::complex<double>> &roots){
    double pivot = (fabs(from)+fabs(to))/4;
    std::uniform_int_distribution<> distrib(from, pivot);
    std::uniform_int_distribution<> shake(shakefrom, shaketo);

    std::complex<double> x1(distrib(gen), distrib(gen));
    std::complex<double> x2(distrib(gen) + shake(gen) + epsilon, distrib(gen) + shake(gen) + epsilon);
    std::complex<double> x3(distrib(gen) + shake(gen) + 2*epsilon, distrib(gen) + shake(gen) + 2*epsilon);

    roots = {x1, x2, x3};
}
*/

int main(){
    // Randomizer
    std::random_device rd;
    std::mt19937 gen(rd());

    std::setprecision(30);
    auto SolverB = Baydoun<double>();
    auto SolverV = Vieta<double>();
    std::vector<std::complex<double>> roots;
    int roots_num = SolverB(2.0, 3.5, 4.0, 5.0, roots);
    std::cout << "Baydoun. Найдено " << roots_num << " корня.\n";
    for (auto &root : roots){
        std::cout << root << "\n";
    }
    roots_num = SolverV(2.0, 3.5, 4.0, 5.0, roots);
    std::cout << "Vieta. Найдено " << roots_num << " корня.\n";
    for (auto &root : roots){
        std::cout << root << "\n";
    }

    double arr[][4] = {
        {12.0, 34.2, 56.99, 7.2},
        {1.0001, 2.03, 3.02, 4.01},
    };
    std::vector<std::vector<std::complex<double>>> roots_mul;
    int *mul_num = SolverB(arr, 2, roots_mul);
    std::cout << "Baydoun. Для 2 полиномов:\n";
    for(int i = 0; i < 2; i++){
        std::cout << "Найдено " << mul_num[i] << " корня.\n";
        std::cout << arr[i][0] << "x^3 + " << arr[i][1] << "x^2 + " << arr[i][2] << "x + " << arr[i][3] << "\n";
        for (auto &root : roots_mul[i]){
           std::cout << root << "\n";
        }
    }

    mul_num = SolverV(arr, 2, roots_mul);
    std::cout << "Vieta. Для 2 полиномов:\n";
    for(int i = 0; i < 2; i++){
        std::cout << "Найдено " << mul_num[i] << " корня.\n";
        std::cout << arr[i][0] << "x^3 + " << arr[i][1] << "x^2 + " << arr[i][2] << "x + " << arr[i][3] << "\n";
        for (auto &root : roots_mul[i]){
           std::cout << root << "\n";
        }
    }
    delete[] mul_num;
    return 0;
}
