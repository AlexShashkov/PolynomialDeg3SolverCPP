#include "Methods.h"
#include <vector>
#include <iomanip>

int main(){
    std::setprecision(30);
    auto Solver = Baydoun<double>();
    std::vector<std::complex<double>> roots;
    int roots_num = Solver(2.0, 3.5, 4.0, 5.0, roots);
    for (auto &root : roots){
        std::cout << root << "\n";
    }
    return 0;
}
