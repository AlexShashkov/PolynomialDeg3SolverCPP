#include "Methods.h"
#include "Generate.h"

int main(){
    std::setprecision(30);
    std::vector<std::complex<double>> gen_roots;
    auto polynomial = generatePolynomial(-1, 1, 0.0001, -0.00001, 0.00001, gen_roots);
    std::cout << "Сгенерированные корни:\n";
    for (auto &root : gen_roots){
       std::cout << root << "\n";
    }
    std::cout << "Коэффициенты полинома X^3, x^2, x, c:\n";
    for(int i = 0; i<4; i++){
        std::cout << polynomial[i] << " ";
    }
    std::cout << "\n";

    Baydoun<double> SolverB;
    auto SolverV = Vieta<double>();
    std::vector<std::complex<double>> roots;
    int roots_num = SolverB(1.0, -1.81214, 0.502814, 0.281353, roots);
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
    delete[] polynomial;
    return 0;
}
