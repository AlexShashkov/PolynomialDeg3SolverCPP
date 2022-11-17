#ifndef METHODS_H
#define METHODS_H
#include <cmath>
#include <vector>
#include <complex>
#include <iostream>
#include <algorithm>
#include <iomanip>

using namespace std::complex_literals;

/*
Имплементация 'Analytical formula for the roots of the general complex cubic polynomial'
Автор: Ibrahim Baydoun
Статья: https://arxiv.org/abs/1512.07585
*/
template <typename number>
class Baydoun
{
	long double _onethree;
	long double _one27;
	long double _sqrt3;
	long double _cbrt4ftwo;


	/*Вычисление вспомогательных степеней коэффициентов полинома.
	@type b: TEMPLATE*
	@param b: Массив, в котором будут храниться коэффициенты x^2.
	@type c: TEMPLATE*
	@param c: Массив, в котором будут храниться коэффициенты x.
	@type d: TEMPLATE*
	@param d: Массив, в котором будут храниться коэффициенты C.
	*/
	void prepare(number *b, number *c, number *d){
		for (int i = 1; i < 6; i++){
		    b[i] = b[i-1]*b[0];
		}
		for (int i = 1; i < 4; i++){
		    c[i] = c[i-1]*c[0];
		}
		for (int i = 1; i < 3; i++){
		    d[i] = d[i-1]*d[0];
		}
	}

	/*Второй шаг из статьи для вычисления корней
	@type b: TEMPLATE*
	@param b: Массив, в котором хранятся коэффициенты x^2.
	@type c: TEMPLATE*
	@param c: Массив, в котором хранятся коэффициенты x.
	@type d: TEMPLATE*
	@param d: Массив, в котором хранятся коэффициенты C.
	@rtype: vector<complex<TEMPLATE>>
	@returns: Корни уравнения.
	*/
	std::vector<std::complex<number>> part2(number *b, number *c, number *d,
			number o, number r,number bthree){
		number onethree = static_cast<number>(_onethree);
		number sqrt3    = static_cast<number>(_sqrt3);
		number cbrt4ftwo= static_cast<number>(_cbrt4ftwo);
		
		number c0d0 = c[0]*d[0];
		number b0c0 = b[0]*c[0];
		number b0c1 = b[0]*c[1];
		number b1c1 = b[1]*c[1]*static_cast<number>(4);
		number t = c[2]*(static_cast<number>(16)*b[5] + static_cast<number>(72)*d[1] + static_cast<number>(264)*b[2]*d[0] + static_cast<number>(66)*b[1]*c[1] - \
		    static_cast<number>(132)*b[0]*c0d0 + static_cast<number>(2)*c[2]) + b[3]*c[0]*(static_cast<number>(12)*d[1] - static_cast<number>(84)*c[2]) - b[1] * \
		    c[1]*d[0]*(static_cast<number>(24)*b[2]+static_cast<number>(291)*d[0]) + d[2]*(static_cast<number>(144)*b0c0 - static_cast<number>(27)*d[0] - static_cast<number>(2)*b[2]);
		number d0 = static_cast<number>(4)*(b[3]*c[1] - b[2]*c0d0) - static_cast<number>(12)*c[0]*d[1] - static_cast<number>(14)*b[1]*c[2] + \
			static_cast<number>(28)*b0c1*d[0] + b[1]*d[1] + c[3];
		std::complex<number> sqrt1;
		if (o > 0)
		    sqrt1 = sqrt(o);
		else
		    sqrt1 = std::complex<number>(0, 1)*static_cast<number>(sqrt(fabs(o)));
		std::complex<number> sqrt2 = std::complex<number>(0, 1)*sqrt3;
		auto sqrt2div3 = sqrt2*onethree;
		auto sqrt2div9 = sqrt2div3*onethree;

		auto bl = (d[0]-b0c0) * sqrt1 * (b1c1 - static_cast<number>(4)*b[0]*c0d0 + static_cast<number>(2)*c[2] + d[1]) +sqrt2div9*t;
		auto bl1 = pow(bl, onethree);
		auto bl2 = pow(bl1, static_cast<number>(2.0));
		auto A1 = (-sqrt2div3)*(static_cast<number>(8)*b[2]*c[0] - static_cast<number>(4)*d[0]*b[1] - static_cast<number>(26)*b0c1 + static_cast<number>(30)*c0d0) + static_cast<number>(2)*c[0]*sqrt1;
		auto A2 = static_cast<number>(8)*(b[4]*c[1] - b[3]*c0d0) - static_cast<number>(40)*b[2]*c[2] + static_cast<number>(2)*b[2]*d[1] +\
		    static_cast<number>(29)*b1c1*d[0] + static_cast<number>(23)*b[0]*c[3] - static_cast<number>(21)*c[2]*d[0] +\
		    static_cast<number>(27)*d[2] - static_cast<number>(99)*b0c0*d[1] -sqrt1*sqrt2 * (static_cast<number>(2)*b1c1 - static_cast<number>(10)*b[0]*c0d0 + c[2] + static_cast<number>(3)*d[1]);
		auto Rbase = sqrt1 * sqrt2div9;
		std::complex<number> R1, R2;
		if (o == 0)
		    if (r > 0){
		R1 = pow(Rbase + r, onethree);
		R2 = -R1;
		    }
		    else{
		R2 = pow(Rbase - r, onethree);
		R1 = -R2;
		    }
		else{
		    R1 = pow(Rbase + r, onethree);
		    R2 = pow(Rbase - r, onethree);
		}
		auto sqrt205=sqrt2*static_cast<number>(0.5);
		std::complex<number> M[2] = {static_cast<number>(0.5) - sqrt205, static_cast<number>(0.5) + sqrt205};
		std::complex<number> M2[2] = {static_cast<number>(-0.5) - sqrt205, static_cast<number>(-0.5) + sqrt205};

		auto arg1_1 = A1*bl1;
		auto arg1_2 = -d0*R1;
		auto arg2_1 = A2*bl2;
		auto arg2_2 = pow(d0, static_cast<number>(2.0))*R2;
		// Вычисляем аргумент комплексного числа
		auto phi1 = std::arg(arg1_1) - std::arg(arg1_2);
		auto phi2 = std::arg(arg2_1) - std::arg(arg2_2);
		auto a1 = (cbrt4ftwo)*(std::cos(phi1)+std::complex<number>(0, std::sin(phi1)));
		auto a2 = (cbrt4ftwo)*(std::cos(phi2)+std::complex<number>(0, std::sin(phi2)));
		auto a1R1 = a1*R1;
		auto a2R2 = a2*R2;
		auto x1 = -a1R1 + a2R2 - bthree;
		auto x2 = M[0]*a1R1 + M2[0]*a2R2 - bthree;
		auto x3 = M[1]*a1R1 + M2[1]*a2R2 - bthree;
		return std::vector<std::complex<number>>{x1, x2, x3};
	}

	/*Начало вычисления.
	@type b: TEMPLATE*
	@param b: Массив, в котором хранятся коэффициенты x^2.
	@type c: TEMPLATE*
	@param c: Массив, в котором хранятся коэффициенты x.
	@type d: TEMPLATE*
	@param d: Массив, в котором хранятся коэффициенты C.
	@type roots: vector<complex<TEMPLATE>>&
	@param root: Вектор, который хранит корни уравнения.
	@rtype: int
	@returns: Количество корней.
	*/
	int solve(number *b, number *c, number *d, std::vector<std::complex<number>> &roots){
        number one27 = static_cast<number>(_one27);
		number bthree = b[0]*static_cast<number>(_onethree);
		number o = -static_cast<number>(4)*(b[2]*d[0] + c[2]) + b[1]*c[1] + static_cast<number>(18)*b[0]*c[0]*d[0] - static_cast<number>(27)*d[1];
		number r = static_cast<number>(2)*b[2]*one27 -static_cast<number>(9)*b[0]*c[0]*one27 + d[0];

		// std::cout << "_SOLVE o r, one27" << o << " " << r << " " << one27 << "\n";

		if(o == 0 && r == 0){
			roots = std::vector<std::complex<number>>{-bthree, -bthree, -bthree};
			return 3;
		}
		else{
			roots = part2(b, c, d, o, r,bthree);
			for (auto &r: roots) {
				// Проверка на нулевые Im
				// std::cout << std::numeric_limits<number>::epsilon() << "\n";
				if(fabs(r.imag()) < fabs(r)*std::numeric_limits<number>::epsilon()) r.imag(0);
			}
			roots.erase(std::unique( roots.begin(), roots.end() ), roots.end());
			return roots.size();
		};
	}

public:
	Baydoun(){
		_onethree = 1.0L/3.0L;
		_one27 = _onethree*_onethree*_onethree;
		_sqrt3 = std::sqrt(3L);
		_cbrt4ftwo = pow(4L,_onethree)*static_cast<number>(0.5);
	}

	/*Функтор для решения уравнения методом Baydoun.
	@type a: TEMPLATE
	@param a: Коэффициент x^3.
	@type b: TEMPLATE
	@param b: Коэффициент x^2.
	@type c: TEMPLATE
	@param c: Коэффициент x.
	@type d: TEMPLATE
	@param d: Коэффициент C.
	@type roots: vector<complex<TEMPLATE>>&
	@param root: Вектор, который хранит корни уравнения.
	@rtype: int
	@returns: Количество корней.
	*/
	int operator()(number a, number b, number c, number d,
			std::vector<std::complex<number>> &roots){
        // x^3, x^2, x, c
		if(a != 0){
		    number _a = 1/a;
		    b *= _a;
		    c *= _a;
		    d *= _a;
		    a = 1;
		}

		number *_b = new number[6];
		number *_c = new number[4];
		number *_d = new number[3];
		_b[0] = b; _c [0]= c; _d[0] = d;
		prepare(_b, _c, _d);

		int result = solve(_b, _c, _d, roots);
		delete[] _b;
		delete[] _c;
		delete[] _d;
		return result;
	}

	/*Функтор для решения уравнения методом Baydoun.
	@type inp: vector<TEMPLATE>&
	@param a: Коэффициенты.
	@type roots: vector<complex<TEMPLATE>>&
	@param root: Вектор, который хранит корни уравнения.
	@rtype: int
	@returns: Количество корней.
	*/
	int operator()(std::vector<number> &inp, std::vector<std::complex<number>> &roots, bool reverse=false){
		if(reverse)
        	return operator()(inp[3], inp[2], inp[1], inp[0], roots);
		else
			return operator()(inp[0], inp[1], inp[2], inp[3], roots);
	}

	/*Функтор для решения уравнений методом Baydoun.
	@type poly: TEMPLATE**
	@param a: Динамический массив размером (coun, 4), где count - количество полиномов.
	@type count: int
	@param count: Количество полиномов.
	@type roots: vector<vector<complex<TEMPLATE>>>&
	@param root: Вектор, который хранит корни уравнения.
	@rtype: int*
	@returns: Количество корней в каждом полиноме.
	*/
	int* operator()(number **poly, int count,
			std::vector<std::vector<std::complex<number>>> &roots){
		// x^3, x^2, x, c
		int *numbers = new int[count];
		for(int i = 0; i < count; i++){
			// std::cout << i << "\n";
			std::vector<std::complex<number>> res;
			numbers[i] = operator()(poly[i][0], poly[i][1], poly[i][2], poly[i][3], res);
			roots.push_back(res);
		}
		return numbers;
	}

	/*Функтор для решения уравнений методом Baydoun.
	@type poly: TEMPLATE[][4]
	@param a: Массив размером (coun, 4), где count - количество полиномов.
	@type count: int
	@param count: Количество полиномов.
	@type roots: vector<vector<complex<TEMPLATE>>>&
	@param root: Вектор, который хранит корни уравнения.
	@rtype: int*
	@returns: Количество корней в каждом полиноме.
	*/
	int* operator()(number poly[][4], int count,
			std::vector<std::vector<std::complex<number>>> &roots){
		// x^3, x^2, x, c
		int *numbers = new int[count];
		for(int i = 0; i < count; i++){
			std::vector<std::complex<number>> res;
			numbers[i] = operator()(poly[i][0], poly[i][1], poly[i][2], poly[i][3], res);
			roots.push_back(res);
		}
		return numbers;
	}
};


/*
Имплементация Тригонометрической формулы Виеты
http://poivs.tsput.ru/ru/Math/Functions/Polynomials/VietaTrigonometricFormula
*/
template <typename number>
class Vieta
{
        long double _pi2div3;
        long double _sqrt3;
        long double _onethree;

	int sign(number val) {
		if( number(0) < val){ return -1; }
		else{ return 1;}
	}

	/*Вырожденный случай.
	@type R: TEMPLATE
	@param R: Вычисленное значение R.
	@type b: TEMPLATE
	@param b: Коэффициент x^2
	@rtype: vector<complex<TEMPLATE>>
	@returns: Вектор, хранящий корни уравнения.
	*/
	std::vector<std::complex<number>> degenerate(number R, number b){
		std::vector<std::complex<number>> roots;
		auto inp2three = b*static_cast<number>(_onethree);
		auto _x = cbrt(R);
		auto x1 = -static_cast<number>(2)*_x-inp2three;
		auto x2 = _x-inp2three;
		std::cout<<"2 действ";
		roots = {x1, x2};
		for (auto &r: roots) {
			// Проверка на нулевые Im
			// std::cout << std::numeric_limits<number>::epsilon() << "\n";
			if(fabs(r.imag()) < fabs(r)*std::numeric_limits<number>::epsilon()) r.imag(0);
		}
		roots.erase(std::unique( roots.begin(), roots.end() ), roots.end());
		return roots;
	}

	/*Действительные корни.
	@type Q: TEMPLATE
	@param Q: Вычисленное значение Q.
	@type Q3: TEMPLATE
	@param Q3: Вычисленное значение Q^3.
	@type R: TEMPLATE
	@param R: Вычисленное значение R.
	@type b: TEMPLATE
	@param b: Коэффициент x^2
	@rtype: vector<complex<TEMPLATE>>
	@returns: Вектор, хранящий корни уравнения.
	*/
	std::vector<std::complex<number>> usual(number Q, number Q3, number R, number b){
		std::vector<std::complex<number>> roots;
		number x1,x2,x3 = 0;
		number pi2div3 = b*static_cast<number>(_pi2div3);
		number inp2three = b*static_cast<number>(_onethree);
		number phi = acos(R/sqrt(Q3))*static_cast<number>(_onethree);
		number sqrtQ = static_cast<number>(2)*sqrt(Q);
		x1 = -sqrtQ*cos(phi)-inp2three;
		x2 = -sqrtQ*cos(phi+pi2div3)-inp2three;
		x3 = -sqrtQ*cos(phi-pi2div3)-inp2three;
		std::cout<<"3 действ";
		roots = {x1, x2, x3};
		for (auto &r: roots) {
			// Проверка на нулевые Im
			// std::cout << std::numeric_limits<number>::epsilon() << "\n";
			if(fabs(r.imag()) < fabs(r)*std::numeric_limits<number>::epsilon()) r.imag(0);
		}
		roots.erase(std::unique( roots.begin(), roots.end() ), roots.end());
		return roots;
	}

	/*Действительные корни.
	@type Q: TEMPLATE
	@param Q: Вычисленное значение Q.
	@type Q3: TEMPLATE
	@param Q3: Вычисленное значение Q^3.
	@type R: TEMPLATE
	@param R: Вычисленное значение R.
	@type b: TEMPLATE
	@param b: Коэффициент x^2
	@rtype: vector<complex<TEMPLATE>>
	@returns: Вектор, хранящий корни уравнения.
	*/
	std::vector<std::complex<number>> complex(number Q, number Q3, number R, number b){
        number onethree = static_cast<number>(_onethree);
        number sqrt3 = static_cast<number>(_sqrt3);
		std::vector<std::complex<number>> roots;
		number x1 = 0;
		std::complex<number>  x2, x3 = 0;
		number inp2three = b*onethree;
		number _phi= 0;
		number T;
		number Tin;
		number sqrtsh;
		number absQ3 = fabs(Q3);
		number sqrtabsQ = sqrt(fabs(Q));
		if(Q > 0){
			number phi = acosh(fabs(R)/sqrt(absQ3))*onethree;
			T = sqrtabsQ*cosh(phi);
			sqrtsh = sqrt3*sqrtabsQ*sinh(phi);
			Tin = T - inp2three;
		}
		else{
			number phi = asinh(fabs(R)/sqrt(absQ3))*onethree;
			T = sqrtabsQ*sinh(phi);
			sqrtsh = sqrt3*sqrtabsQ*cosh(phi);
			Tin = T - inp2three;
		}
		x1 = static_cast<number>(-2)*T-inp2three;
		x2 = std::complex<number>(Tin,sqrtsh);
		x3 = std::complex<number>(Tin,-sqrtsh);
		roots = {x1, x2, x3};
		for (auto &r: roots) {
			// Проверка на нулевые Im
			// std::cout << std::numeric_limits<number>::epsilon() << "\n";
			if(fabs(r.imag()) < fabs(r)*std::numeric_limits<number>::epsilon()) r.imag(0);
		}
		roots.erase(std::unique( roots.begin(), roots.end() ), roots.end());
		return roots;
	}
public:
	Vieta(){
		_pi2div3 = 2L*M_PI/3L;
		_sqrt3 = sqrt(3L);
		_onethree = 1.0L/3.0L;
	}

        /*Функтор для решения уравнения методом Baydoun.
        @type a: TEMPLATE
        @param a: Коэффициент x^3.
        @type b: TEMPLATE
        @param b: Коэффициент x^2.
        @type c: TEMPLATE
        @param c: Коэффициент x.
        @type d: TEMPLATE
        @param d: Коэффициент C.
        @type roots: vector<complex<TEMPLATE>>&
        @param root: Вектор, который хранит корни уравнения.
        @rtype: int
        @returns: Количество корней.
        */
	int operator()(number a, number b, number c, number d,
			std::vector<std::complex<number>> &roots){
		// x^3, x^2, x, c
		number onethree = static_cast<number>(_onethree);
		if(a != 0){
			number _a = 1/a;
			b *= _a;
			c *= _a;
			d *= _a;
			a = 1;
		}

		auto Q = b*b*onethree*onethree -c*onethree;
		// std::cout << "b " << b << " b^2 " << b*b << " c " << c << "\n";
		// std::cout << "onethree " << onethree << "\n";
		// std::cout << Q << ":)\n";
		if(Q == 0){
			number rs = -b/3;
			roots = {rs, rs, rs};
			return 3;
		}
		else{
			number R = pow(b, 3)*onethree*onethree*onethree-c*b/6+d*0.5;
			auto R2 = R*R;
			auto Q3 = Q*Q*Q;
			auto S = Q3-R2;
			
			if(S==0){
				roots = degenerate(R, b);
			}
			else if(S > 0){
				roots = usual(Q, Q3, R, b);
			}
			else{
				roots = complex(Q, Q3, R, b);
			}
			return roots.size();
		}
	}

	/*Функтор для решения уравнения методом Виета.
	@type inp: vector<TEMPLATE>&
	@param a: Коэффициенты.
	@type roots: vector<complex<TEMPLATE>>&
	@param root: Вектор, который хранит корни уравнения.
	@rtype: int
	@returns: Количество корней.
	*/
	int operator()(std::vector<number> &inp, std::vector<std::complex<number>> &roots, bool reverse=false){
		if(reverse)
        	return operator()(inp[3], inp[2], inp[1], inp[0], roots);
		else
			return operator()(inp[0], inp[1], inp[2], inp[3], roots);
	}

        /*Функтор для решения уравнений методом Виета.
        @type poly: TEMPLATE**
        @param a: Динамический массив размером (coun, 4), где count - количество полиномов.
        @type count: int
        @param count: Количество полиномов.
        @type roots: vector<vector<complex<TEMPLATE>>>&
        @param root: Вектор, который хранит корни уравнения.
        @rtype: int*
        @returns: Количество корней в каждом полиноме.
        */
	int* operator()(number **poly, int count,
			std::vector<std::vector<std::complex<number>>> &roots){
		// x^3, x^2, x, c
		int *numbers = new int[count];

		std::vector<std::complex<number>> res;
		for(int i = 0; i < count; i++){
			// std::cout << i << "\n";
			numbers[i] = operator()(poly[i][0], poly[i][1], poly[i][2], poly[i][3], res);
			roots.push_back(res);
			res.clear();
		}
		return numbers;
	}

        /*Функтор для решения уравнений методом Виета.
        @type poly: TEMPLATE[][4]
        @param a: Массив размером (coun, 4), где count - количество полиномов.
        @type count: int
        @param count: Количество полиномов.
        @type roots: vector<vector<complex<TEMPLATE>>>&
        @param root: Вектор, который хранит корни уравнения.
        @rtype: int*
        @returns: Количество корней в каждом полиноме.
        */
	int* operator()(number poly[][4], int count,
			std::vector<std::vector<std::complex<number>>> &roots){
		// x^3, x^2, x, c
		int *numbers = new int[count];
		std::vector<std::complex<number>> res;
		for(int i = 0; i < count; i++){
			numbers[i] = operator()(poly[i][0], poly[i][1], poly[i][2], poly[i][3], res);
			roots.push_back(res);
			res.clear();
		}
		
		return numbers;
	}
};
#endif
