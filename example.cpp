#include "Methods.h"
#include "Generate.h"

#define number double

int main(){
    std::setprecision(30);
    std::vector<std::complex<number>> gen_roots;
    auto polynomial = generatePolynomial<number>(-1, 1, 0.0001, -0.00001, 0.00001, gen_roots);
    std::cout << "Сгенерированные корни:\n";
    for (auto &root : gen_roots){
       std::cout << root << "\n";
    }
    std::cout << "Коэффициенты полинома X^3, x^2, x, c:\n";
    for(int i = 0; i<4; i++){
        std::cout << polynomial[i] << " ";
    }
    std::cout << "\n";

    Baydoun<number> SolverB;
    auto SolverV = Vieta<number>();
    std::vector<std::complex<number>> roots;
    int roots_num = SolverB(polynomial[0],polynomial[1],polynomial[2],polynomial[3], roots);
    std::cout << "Baydoun. Найдено " << roots_num << " корня.\n";
    for (auto &root : roots){
        std::cout << root << "\n";
    }
    roots_num = SolverV(polynomial[0],polynomial[1],polynomial[2],polynomial[3], roots);
    std::cout << "Vieta. Найдено " << roots_num << " корня.\n";
    for (auto &root : roots){
        std::cout << root << "\n";
    }
    number arr[][4] = {
        {12.0, 34.2, 56.99, 7.2},
        {1.0, 1.0, 1.0, 1.0},
    };
    std::vector<std::vector<std::complex<number>>> roots_mul;
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
