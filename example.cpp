#include "Methods.h"
#include "Generate.h"

#define number float

int main(){
    std::setprecision(30);

    Baydoun<number> SolverB;
    auto SolverV = Vieta<number>();
    std::vector<std::complex<number>> roots;
    // x^3 + (0.485677)x^2 + (0.0786274)x + (0.00424306) = 0
    int roots_num = SolverB(1.0f, 0.485677f, 0.0786274f, 0.00424306f, roots);
    std::cout << "Baydoun. Найдено " << roots_num << " корня.\n";
    for (auto &root : roots){
        std::cout << root << "\n";
    }
    
    roots_num = SolverV(1.0f, 0.485677f, 0.0786274f, 0.00424306f, roots);
    std::cout << "Vieta. Найдено " << roots_num << " корня.\n";
    for (auto &root : roots){
        std::cout << root << "\n";
    }
    return 0;
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
    return 0;
}
