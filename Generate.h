#include <iostream>
#include <complex>
#include <vector>
#include <random>
#include <cmath>

/*Генерация корней
@type from: double
@param from: от
@type to: double
@param to: до
@type epsilon: double
@param epsilon: расстояние между корнями
@type shakefrom: double
@param shakefrom: смещение корня от
@type shaketo: double
@param shaketo: смещение корня до
@type gen: mt19937
@param gen: Генератор чисел
@rtype: vector<complex<double>>
@returns: Сгенерированные корни
*/
std::vector<std::complex<double>> generateEpsilon(double from, double to, double epsilon, double shakefrom, double shaketo, std::mt19937 gen){
    std::uniform_real_distribution<> distrib(from, to - epsilon);
    std::uniform_real_distribution<> shake(shakefrom, shaketo);

    auto init = distrib(gen);
    std::complex<double> x1(init, 0);
    auto a = init + shake(gen) + epsilon;
    auto b = init + shake(gen) + epsilon;
    std::complex<double> x2(a, b);
    std::complex<double> x3(a, -b);

    return std::vector<std::complex<double>>{x1, x2, x3};
}

/*Генерация коэффициента
@type k: int
@param root: Позиция коэффициента от 1 до 3 (x^2 до C)
@type roots: vector<complex<TEMPLATE>>&
@param root: Вектор, который хранит корни уравнения.
@rtype: double
@returns: Сгенерированный коэффициент
*/
double coeff(int k, std::vector<std::complex<double>>& roots)
{
    int size = roots.size();
    std::vector<int> loop_counter(k+1); // nested loops

    for(int i=0; i<k+1; ++i)
        loop_counter[i]=0;


    std::vector<int> max(k+1);
    for(int i=0 ; i<k ; ++i)
        max[i] = size-k+i+1;

    max[k] = (int)INFINITY;
    int p1 = 0;
    std::complex<double> sum = 0;

    int counter;
    while(loop_counter[k]==0)
    {
        // variable nested loops
        counter = 0;
        for(int i = 1 ; i < k; ++i)
        {
            if(loop_counter[i-1] < loop_counter[i])
            ++counter;
        }

        if(counter == k - 1) // true if i_1 < i_2 < i_3 ....
        {
            std::complex<double> prod(1);
            for(int i = 0 ; i < k ; ++i)
            prod *= roots[loop_counter[i]];   // taking products

            sum += prod;  // increament
        }

        ++loop_counter[0];
        while(loop_counter[p1]==max[p1]){
            loop_counter[p1]=0;
            loop_counter[++p1]++;
            if(loop_counter[p1]!=max[p1])
                p1=0;
        }
    }
    return pow(-1.0,k)*sum.real();
}

/*Генерация коэффициентов и корней
@type from: double
@param from: от
@type to: double
@param to: до
@type epsilon: double
@param epsilon: расстояние между корнями
@type shakefrom: double
@param shakefrom: смещение корня от
@type shaketo: double
@param shaketo: смещение корня до
@type roots: vector<complex<TEMPLATE>>&
@param root: Вектор, который хранит корни уравнения.
@rtype: double*
@returns: Сгенерированные коэффициенты
*/
double* generatePolynomial(double from, double to, double epsilon, double shakefrom, double shaketo, std::vector<std::complex<double>>& roots){
    // Randomizer
    std::random_device rd;
    std::mt19937 gen(rd());
    roots = generateEpsilon(from, to, epsilon, shakefrom, shaketo, gen);
    double *coeffs = new double[4];
    coeffs[0] = 1;
    for(int i = 1; i<4; i++){
        coeffs[i] = coeff(i, roots);
    }
    return coeffs;
}

