#include <iostream>
#include <complex>
#include <vector>
#include <random>
#include <cmath>

/*Генерация корней
@type from: number
@param from: от
@type to: number
@param to: до
@type epsilon: number
@param epsilon: расстояние между корнями
@type shakefrom: number
@param shakefrom: смещение корня от
@type shaketo: number
@param shaketo: смещение корня до
@type gen: mt19937
@param gen: Генератор чисел
@rtype: vector<complex<number>>
@returns: Сгенерированные корни
*/
template<typename number>
std::vector<std::complex<number>> generateEpsilon(number from, number to, number epsilon, number shakefrom, number shaketo, std::mt19937 gen){
    std::uniform_real_distribution<> distrib(from, to - epsilon);
    std::uniform_real_distribution<> shake(shakefrom, shaketo);

    auto init = distrib(gen);
    std::complex<number> x1(init, 0);
    auto a = init + shake(gen) + epsilon;
    auto b = init + shake(gen) + epsilon;
    std::complex<number> x2(a, b);
    std::complex<number> x3(a, -b);

    return std::vector<std::complex<number>>{x1, x2, x3};
}

/*Генерация коэффициента
@type k: int
@param root: Позиция коэффициента от 1 до 3 (x^2 до C)
@type roots: vector<complex<TEMPLATE>>&
@param root: Вектор, который хранит корни уравнения.
@rtype: number
@returns: Сгенерированный коэффициент
*/
template<typename number>
number coeff(int k, std::vector<std::complex<number>>& roots)
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
    std::complex<number> sum = 0;

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
            std::complex<number> prod(1);
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
@type from: number
@param from: от
@type to: number
@param to: до
@type epsilon: number
@param epsilon: расстояние между корнями
@type shakefrom: number
@param shakefrom: смещение корня от
@type shaketo: number
@param shaketo: смещение корня до
@type roots: vector<complex<TEMPLATE>>&
@param root: Вектор, который хранит корни уравнения.
@rtype: number*
@returns: Сгенерированные коэффициенты
*/
template<typename number>
number* generatePolynomial(number from, number to, number epsilon, number shakefrom, number shaketo, std::vector<std::complex<number>>& roots){
    // Randomizer
    std::random_device rd;
    std::mt19937 gen(rd());
    roots = generateEpsilon(from, to, epsilon, shakefrom, shaketo, gen);
    number *coeffs = new number[4];
    coeffs[0] = 1;
    for(int i = 1; i<4; i++){
        coeffs[i] = coeff(i, roots);
    }
    return coeffs;
}

