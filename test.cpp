#include "Methods.h"
#include "Generate.h"
#include <fstream>

bool compareComplexVec(std::vector<std::complex<double>> a, std::vector<std::complex<double>> b, double allowed){
    for (int i = 0; i < a.size(); i++){
        if (std::fabs(a[i] - b[i]) >= allowed) return false;
    }
    return true;
}

int main(int argc,char* argv[]){
    std::setprecision(30);
    int count;
    double allowed = 1e-2;
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
    else if (argc == 3){
        count = std::atoi(argv[1]);
        allowed  = std::atof(argv[2]);
    }
    else if (argc == 6){
        count = std::atoi(argv[1]);
        allowed  = std::atof(argv[2]);
        from  = std::atof(argv[2]);
        to    = std::atof(argv[3]);
        step  = std::atof(argv[4]);
    }
    else if (argc == 8){
        count        = std::atoi(argv[1]);
        allowed      = std::atof(argv[2]);
        from         = std::atof(argv[3]);
        to           = std::atof(argv[4]);
        step         = std::atof(argv[5]);
        epsilon_from = std::atof(argv[6]);
        epsilon_to   = std::atof(argv[7]);
    }
    std::cout << "Launching " << count << " " << allowed << " " << from << " " << to << " " 
        << step << " " << epsilon_from << " " << epsilon_to << "\n";
    std::cout << "Генерирую полиномы.\n";
    
    double** polynomials = new double*[count];
    std::vector<std::vector<std::complex<double>>> roots;
    for(int i = 0; i < count; i++){
        std::vector<std::complex<double>> gen;
        polynomials[i] = generatePolynomial(from, to, step, epsilon_from, epsilon_to, gen);
        roots.push_back(gen);
    }

    // Сортировка корней
    auto csort = [](std::complex<double> a, std::complex<double> b){
        if (imag(a) == imag(b))
            return real(a) < real(b);
        return imag(a) < imag(b);
    };

    for(int i = 0; i < count; i++){
        sort(
            roots[i].begin(),
            roots[i].end(),
            csort
        );
    }

    auto SolverB = Baydoun<double>();
    auto SolverV = Vieta<double>();

    std::vector<std::vector<std::complex<double>>> rootsB;
    std::vector<std::vector<std::complex<double>>> rootsV;
    std::cout << "Решаю с помощью Baydoun.\n";
    auto countRB = SolverB(polynomials, count, rootsB);
    std::cout << "Решаю с помощью Виета.\n";
    auto countRV = SolverV(polynomials, count, rootsV);

    for(int i = 0; i < count; i++){
        sort(
            rootsB[i].begin(),
            rootsB[i].end(),
            csort
        );
        sort(
            rootsV[i].begin(),
            rootsV[i].end(),
            csort
        );
    }

    std::cout << "\nНачинаю запись полиномов и результатов в файл\n";

    std::ofstream file;
    file.open ("polynomials.txt");
    file << count << " " << from << " " << to << " " << step << " " << epsilon_from << " " << epsilon_to << "\n";
    for(int i = 0; i < count; i++){
        for(int j = 0; j < 4; j++)
          file << polynomials[i][j] << " ";
        file << "; ";
        for (auto &root : roots[i]){
           file << root << " ";
        }
        file << "\n";
    }
    file.close();

    file.open ("vietaresults.txt");
    file << "Vieta results\n";
    for(int i = 0; i < count; i++){
        for (auto &root : rootsV[i]){
           file << root << " ";
        }
        file << compareComplexVec(roots[i], rootsV[i], allowed) << "\n";
    }
    file.close();

    file.open ("baydounresults.txt");
    file << "Baydoun results\n";
    for(int i = 0; i < count; i++){
        for (auto &root : rootsB[i]){
           file << root << " ";
        }
        file << compareComplexVec(roots[i], rootsB[i], allowed) << "\n";
    }
    file.close();
    std::cout << "Результаты сохранены.\n";
    for(int i = 0; i < count; ++i) {
        delete [] polynomials[i];
    }
    delete [] polynomials;
    delete [] countRB;
    delete [] countRV;
}
