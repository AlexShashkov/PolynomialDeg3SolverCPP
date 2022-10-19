#include "Methods.h"
#include <vector>

int main(){
    auto Solver = Baydoun<long double>();
    std::vector<std::complex<long double>> roots;
    int roots_num = Solver(2.0L, 3.5L, 4.0L, 5.0L, roots);
    for (auto &root : roots){
        std::cout << root << "\n";
    }
    return 0;
}
