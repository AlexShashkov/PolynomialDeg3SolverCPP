#include "Methods.h"
#include "Generate.h"

int main(int argc,char* argv[]){
    std::setprecision(30);
    int count;
    double from = -1;
    double to = 1;
    double step = 0.1;
    double epsilon_from = 0.01;
    double epsilon_to = 0.01;

    if (argc == 1){
        std::cout << "Введите количество полиномов > ";
        std::cin >> count;
    }
    else if (argc == 2){
        count = std::atoi(argv[1]);
    }
    else if (argc == 4){
        count = std::atoi(argv[1]);
        from  = std::atof(argv[2]);
        to    = std::atof(argv[3]);
    }
    else if (argc == 5){
        count = std::atoi(argv[1]);
        from  = std::atof(argv[2]);
        to    = std::atof(argv[3]);
        step  = std::atof(argv[4]);
    }
    else if (argc == 7){
        count        = std::atoi(argv[1]);
        from         = std::atof(argv[2]);
        to           = std::atof(argv[3]);
        step         = std::atof(argv[4]);
        epsilon_from = std::atof(argv[5]);
        epsilon_to   = std::atof(argv[6]);
    }
    std::cout << "Генерирую полиномы.\n";
    
    double** polynomials = new double*[count];
    std::vector<std::vector<std::complex<double>>> roots;
    for(int i = 0; i < count; i++){
        std::vector<std::complex<double>> gen;
        polynomials[i] = generatePolynomial(from, to, step, epsilon_from, epsilon_to, gen);
        roots.push_back(gen);
    }

    std::cout << "Сгенерированные полиномы.\n";
    std::cout << "X^3, x^2, x, C\n";
    for(int i = 0; i < count; i++){
        for(int j = 0; j < 4; j++)
          std::cout << polynomials[i][j] << " ";
        std::cout << "\n";
    }

    auto SolverB = Baydoun<double>();
    auto SolverV = Vieta<double>();

    std::vector<std::vector<std::complex<double>>> rootsB;
    std::vector<std::vector<std::complex<double>>> rootsV;
    std::cout << "Решаю с помощью Baydoun.\n";
    auto countRB = SolverB(polynomials, count, rootsB);
    std::cout << "Решаю с помощью Виета.\n";
    auto countRV = SolverV(polynomials, count, rootsV);

    std::cout << "Результаты Baydoun:\n";
    for(int i = 0; i < count; i++){
        std::cout << "Найдено " << countRB[i] << " корня.\n";
        std::cout << polynomials[i][0] << "x^3 + " << polynomials[i][1] << "x^2 + " << polynomials[i][2] << "x + " << polynomials[i][3] << "\n";
        std::cout << "Ответ: ";
        for (auto &root : rootsB[i]){
           std::cout << root << " ";
        }
        std::cout << "\nОжидалось: ";
        for (auto &root : roots[i]){
           std::cout << root << " ";
        }
    }

    std::cout << "\n\nРезультаты Виета:\n";
    for(int i = 0; i < count; i++){
        std::cout << "Найдено " << countRV[i] << " корня.\n";
        std::cout << polynomials[i][0] << "x^3 + " << polynomials[i][1] << "x^2 + " << polynomials[i][2] << "x + " << polynomials[i][3] << "\n";
        std::cout << "Ответ: ";
        for (auto &root : rootsV[i]){
           std::cout << root << " ";
        }
        std::cout << "\nОжидалось: ";
        for (auto &root : roots[i]){
           std::cout << root << " ";
        }
    }
    std::cout << "\n";
    for(int i = 0; i < count; ++i) {
        delete [] polynomials[i];
    }
    delete [] polynomials;
    delete [] countRB;
    delete [] countRV;
}
