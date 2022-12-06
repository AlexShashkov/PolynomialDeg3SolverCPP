#ifndef METHODS_H
#define METHODS_H

#include <cmath>
#include <numbers>
#include <vector>
#include <complex>
#include <iostream>
#include <algorithm>
#include <iomanip>

using namespace std::complex_literals;

namespace implementations{
	/*Аргумент комплексного числа
	@type inp: complex<TEMPLATE>
	@param inp: Комплексное число.
	@rtype: TEMPLATE
	@returns: Аргумент комплексного числа.
	*/
	template <typename number>
	number argp(std::complex<number> inp){
		const number _PI = std::numbers::pi_v<number>;
		number x = std::real(inp);
		number y = std::imag(inp);
		if(x > 0) return std::arg(inp);
		else{
			number _pi = y < 0 ? -_PI: _PI;
			if(x == 0) return -_pi;
			return std::arg(inp) + _pi;
		}
	}

	/*FMS
	@type a: TEMPLATE
	@type b: TEMPLATE
	@type c: TEMPLATE
	@type d: TEMPLATE
	@rtype:  TEMPLATE
	*/
	template <typename number>
	number fms(number a, number b, number c, number d) {
		auto tmp = d * c;
		return std::fma(a, b, -tmp);
	}

	/*FMS
	@type a: complex<TEMPLATE>
	@type b: complex<TEMPLATE>
	@type c: complex<TEMPLATE>
	@type d: complex<TEMPLATE>
	@rtype:  complex<TEMPLATE>
	*/
	template <typename number>
	std::complex<number> fms(std::complex<number> a, std::complex<number> b, std::complex<number> c,std::complex<number> d) {
		return {fms(a.real(),b.real(),c.real(),d.real())+fms(c.imag(),d.imag(),a.imag(),b.imag()),
		fms(a.real(),b.imag(),c.real(),d.imag())+fms(b.real(),a.imag(),d.real(),c.imag())};
	}

	/*fused multiply-add for complex numbers
	@type a: complex<TEMPLATE>
	@type b: complex<TEMPLATE>
	@type c: complex<TEMPLATE>
	@returns: a*b+c
	@rtype:  complex<TEMPLATE>
	*/
	template <typename number>
	std::complex<number> fma(std::complex<number> a, std::complex<number> b, std::complex<number> c) {
		return {fms(a.real(), b.real(), a.imag(), b.imag()) + c.real(),
			fma(a.real(), b.imag(), fma(a.imag(), b.real(), c.imag()))
		};
	}

	/*
	Имплементация 'Analytical formula for the roots of the general complex cubic polynomial'
	Автор: Ibrahim Baydoun
	Статья: https://arxiv.org/abs/1512.07585
	*/
	template <typename number>
	class Baydoun
	{
		const number PI = std::numbers::pi_v<number>;
		number onethree;
		number one27;
		number sqrt3;
		number cbrt4ftwo;
		number THRESHOLD;

		/*Вычисление вспомогательных степеней коэффициентов полинома.
		@type b: TEMPLATE*
		@param b: Массив, в котором будут храниться коэффициенты x^2.
		@type c: TEMPLATE*
		@param c: Массив, в котором будут храниться коэффициенты x.
		@type d: TEMPLATE*
		@param d: Массив, в котором будут храниться коэффициенты C.
		*/
		void prepare(number *b, number *c, number *d){
			number b0 = b[0];
			number c0 = c[0];
			number d0 = d[0];
			for (int i = 1; i < 6; i++){
				b[i] = b[i-1]*b0;
			}
			for (int i = 1; i < 4; i++){
				c[i] = c[i-1]*c0;
			}
			for (int i = 1; i < 3; i++){
				d[i] = d[i-1]*d0;
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
				number o, number r, number bthree, number tmp1){
			// Самые часто вызываемые переменные

			number b0 = b[0];
			number b1 = b[1];
			number b2 = b[2];

			number c0 = c[0];
			number c1 = c[1];
			number c2 = c[2];

			number d0 = d[0];
			number d1 = d[1];
			
			number tmp = std::fma(-b0,c0,d0);

			number b0c0 = -tmp + d0;
			number b0c1 = b0*c1;
			number b1c1 = b1*c1;

			number t = fms(c2, fms(static_cast<number>(2), std::fma(static_cast<number>(12)*d0, 
				fms(static_cast<number>(3), d0, static_cast<number>(-11), b2), std::fma(static_cast<number>(8), b[5], c2)),
				static_cast<number>(66)*b0c0, (tmp+d0)), static_cast<number>(12)*b[3]*c0, std::fma(static_cast<number>(7),c2,-d1)) + \
				fms(d[2], fms(static_cast<number>(2)*b0, std::fma(static_cast<number>(72), c0, -b1),static_cast<number>(27),d0),
				b1c1, d0*static_cast<number>(3)*fms(static_cast<number>(8), b2, static_cast<number>(-97),d0));
			number partiond0 = std::fma(tmp,fms(static_cast<number>(14)*b0,c1,static_cast<number>(4)*b2,c0),fms(static_cast<number>(14)*b0c1,d0,static_cast<number>(12)*c0,d1))+ std::fma(b1,d1,c[3]);
			std::complex<number> sqrt1;
			if (o > THRESHOLD)
				sqrt1 = sqrt(o);
			else
				sqrt1 = std::complex<number>(0, 1)*static_cast<number>(sqrt(fabs(o)));
			std::complex<number> sqrt2 = std::complex<number>(0, 1)*sqrt3;
			std::complex<number> sqrt2div3 = sqrt2*onethree;
			std::complex<number> sqrt2div9 = sqrt2div3*onethree;

			std::complex<number> bl = tmp * sqrt1 * std::fma(static_cast<number>(4)*b0c0,-tmp,static_cast<number>(2)*c2+d1) + sqrt2div9*t;
			std::complex<number> bl1 = pow(bl, onethree);
			std::complex<number> bl2 = pow(bl1, static_cast<number>(2.0));
			std::complex<number> A1 = (-sqrt2div3)*static_cast<number>(2)*std::fma(b1,std::fma(static_cast<number>(2),d0,tmp1),fms(static_cast<number>(15), c0*d0, static_cast<number>(13),b0c1)) + static_cast<number>(2)*c0*sqrt1;
			std::complex<number> A2 = fms(static_cast<number>(2)*b1*d0, fms(b0, d0, static_cast<number>(-58), c1), static_cast<number>(8)*b2,
				fms(b0c0, tmp, static_cast<number>(-5), c2)) + \
				fms(b0c0, fms(static_cast<number>(23), c2, static_cast<number>(99), d1), static_cast<number>(-3)*d0, fms(static_cast<number>(9), d1, 
				static_cast<number>(7),c2)) - sqrt1*sqrt2 *std::fma(static_cast<number>(2)*b0c0,fms(static_cast<number>(4),b0c0,static_cast<number>(5),d0),
				std::fma(static_cast<number>(3),d1,c2));
			std::complex<number>  Rbase = sqrt1 * sqrt2div9;
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
			std::complex<number>  sqrt205=sqrt2*static_cast<number>(0.5);
			std::complex<number> M[2] = {static_cast<number>(0.5) - sqrt205, static_cast<number>(0.5) + sqrt205};

			std::complex<number>  arg1_1 = A1*bl1;
			std::complex<number>  arg1_2 = -partiond0*R1;
			std::complex<number>  arg2_1 = A2*bl2;
			std::complex<number>  arg2_2 = static_cast<number>(pow(partiond0, static_cast<number>(2)))*R2;
			// Вычисляем аргумент комплексного числа
			number phi1 = argp(arg1_1) - argp(arg1_2);
			number phi2 = argp(arg2_1) - argp(arg2_2); 
			std::complex<number> a1 = (cbrt4ftwo)*(std::cos(phi1)+std::complex<number>(0, std::sin(phi1)));
			std::complex<number> a2 = (cbrt4ftwo)*(std::cos(phi2)+std::complex<number>(0, std::sin(phi2)));
			std::complex<number> x1 = fms(a2,R2,a1,R1) - bthree;
			std::complex<number> x2 = fms(M[0]*a1,R1,M[1]*a2,R2) - bthree;
			std::complex<number> x3 = fms(M[1]*a1,R1,M[0]*a2,R2) - bthree;
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
			number d0 = d[0];
			number b0 = b[0];
			number c0 = c[0];

			number tmp1 = std::fma(static_cast<number>(2)*b0,c0,static_cast<number>(-3)*d0);

			number o = std::fma(-b[1],std::fma(static_cast<number>(4)*b0,d0,-c[1]),fms(static_cast<number>(9)*d0,tmp1, static_cast<number>(4),c[2]));
			number r = std::fma(one27,fms(static_cast<number>(2), b[2], static_cast<number>(9), b0*c0),d0);
			number bthree = b0*onethree;

			if(o == 0 && r == 0){
				roots = std::vector<std::complex<number>>{-bthree, -bthree, -bthree};
				return 3;
			}
			else{
				roots = part2(b, c, d, o, r,bthree,tmp1);
				return roots.size();
			};
		}

	public:
		Baydoun(number trsh = 1e-5){
			long double _onethree = 1.0L/3.0L;
			long double _one27 = _onethree*_onethree*_onethree;
			long double _sqrt3 = std::sqrt(3L);
			long double _cbrt4ftwo = pow(4L,_onethree)*static_cast<number>(0.5);

			onethree = static_cast<number>(_onethree);
			one27 = static_cast<number>(_one27);
			sqrt3 = static_cast<number>(_sqrt3);
			cbrt4ftwo = static_cast<number>(_cbrt4ftwo);
			THRESHOLD=trsh;
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
			if(a != 0 && std::isfinite(a)){
				b /= a;
				c /= a;
				d /= a;
				a = 1;
			}
			else
				throw std::invalid_argument("Коэффициент при x^3 равен нулю или б/м.");

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
		const number PI = std::numbers::pi_v<number>;
		number pi2div3;
		number sqrt3;
		number onethree;

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
		std::vector<std::complex<number>> degenerate(number R, number b, number inp2three){
			std::vector<std::complex<number>> roots;
			number _x = cbrt(R);
			auto x1 = -std::fma(static_cast<number>(2),_x,inp2three);
			auto x2 = _x-inp2three;
			roots = {x1, x2};
			for (auto &r: roots) {
				// Проверка на нулевые Im
				if(fabs(r.imag()) < fabs(r)*std::numeric_limits<number>::epsilon()) r.imag(0);
			}
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
		std::vector<std::complex<number>> usual(number Q, number Q3, number R, number b, number inp2three){
			std::vector<std::complex<number>> roots;
			number x1,x2,x3 = 0;
			number phi = acos(R/sqrt(Q3))*onethree;
			number sqrtQ = static_cast<number>(-2)*sqrt(Q);
			x1 = std::fma(sqrtQ,cos(phi),-inp2three);
			x2 = std::fma(sqrtQ,cos(phi+pi2div3),-inp2three);
			x3 = std::fma(sqrtQ,cos(phi-pi2div3),-inp2three);
			roots = {x1, x2, x3};
			for (auto &r: roots) {
				if(fabs(r.imag()) < fabs(r)*std::numeric_limits<number>::epsilon()) r.imag(0);
			}
			return roots;
		}

		/*Комплексные корни.
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
		std::vector<std::complex<number>> complex(number Q, number Q3, number R, number b, number inp2three){
			std::vector<std::complex<number>> roots;
			number x1 = 0;
			std::complex<number>  x2, x3 = 0;
			number _phi= 0;
			number T;
			number Tin;
			number sqrtsh;
			number sqrtabsQ3 = sqrt(fabs(Q3));
			number sqrtabsQ = sqrt(fabs(Q));

			// Q может быть не б/м, но корень может быть еще меньше
			if(sqrtabsQ3 == 0 || !std::isfinite(sqrtabsQ3)){
				return {-inp2three, -inp2three, -inp2three};
			}
			if(Q > 0 && std::isfinite(Q)){
				number phi = acosh(fabs(R)/sqrtabsQ3)*onethree;
				T = sqrtabsQ*cosh(phi);
				sqrtsh = sqrt3*sqrtabsQ*sinh(phi);
			}
			else{
				number phi = asinh(fabs(R)/sqrtabsQ3)*onethree;
				T = sqrtabsQ*sinh(phi);
				sqrtsh = sqrt3*sqrtabsQ*cosh(phi);
			}
			Tin = T - inp2three;
			x1 = -std::fma(static_cast<number>(2),T,inp2three);
			x2 = std::complex<number>(Tin,sqrtsh);
			x3 = std::complex<number>(Tin,-sqrtsh);
			roots = {x1, x2, x3};
			for (auto &r: roots) {
				if(fabs(r.imag()) < fabs(r)*std::numeric_limits<number>::epsilon()) r.imag(0);
			}
			return roots;
		}
	public:
		Vieta(){
			long double _pi2div3 = 2L*PI/3L;
			long double _sqrt3 = sqrt(3L);
			long double _onethree = 1.0L/3.0L;

			pi2div3 = static_cast<number>(_pi2div3);
			sqrt3 = static_cast<number>(_sqrt3);
			onethree = static_cast<number>(_onethree);
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
			if(a != 0 && std::isfinite(a)){
				b /= a;
				c /= a;
				d /= a;
				a = 1;
			}
			else throw std::invalid_argument("Коэффициент при x^3 равен нулю или б/м.");

			number b_onethree = b*onethree;
			auto Q = fms(b_onethree,b_onethree,c,onethree);
			if(Q == 0 || !std::isfinite(Q)){
				number rs = -b_onethree;
				roots = {rs, rs, rs};
				return 3;
			}
			else{
				number R = std::fma(static_cast<number>(0.5), std::fma(-c, b_onethree, d), pow(b_onethree, 3));
				auto Q3 = Q*Q*Q;
				auto S = std::fma(-R,R,Q3);
				if(S==0 || !std::isfinite(S)){
					roots = degenerate(R, b, b_onethree);
				}
				else if(S > 0){
					roots = usual(Q, Q3, R, b, b_onethree);
				}
				else{
					roots = complex(Q, Q3, R, b, b_onethree);
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
}
#endif