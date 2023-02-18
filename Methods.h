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
	using std::complex;
	using std::vector;
	using std::fma;

	/** \brief Check if at least one number is not finite
	*/
	template<typename ... T>
	inline bool anynotfinite(T && ... t)
	{
		return ((!std::isfinite(t)) || ...);
	}

	/** \brief Sign
	*/
	template <typename number>
	inline int sign(number val) {
		return (number(0) < val) - (val < number(0));
	}

	/** \brief Fused multiply-substract
	*/
	template <typename number>
	inline number fms(number a, number b, number c, number d) {
		auto tmp = -d * c;
		return fma(a, b, tmp) + fma(d, c, tmp);
	}

	/** \brief Fused multiply-substract for complex numbers
	*/
	template <typename number>
	inline complex<number> fms(std::complex<number> a, std::complex<number> b, std::complex<number> c,std::complex<number> d) {
		return {
			fms(a.real(),b.real(),c.real(),d.real())+fms(c.imag(),d.imag(),a.imag(),b.imag()),
			fms(a.real(),b.imag(),c.real(),d.imag())+fms(b.real(),a.imag(),d.real(),c.imag())
		};
	}

	/** \brief Fused multiply-substract for complex numbers
	 * 	\return a*b+c
	*/
	template<typename number>
	inline complex<number> fma(complex<number> a, complex<number> b, complex<number> c) {
    return {
			fms(a.real(), b.real(), a.imag(), b.imag()) + c.real(),
            fma(a.real(), b.imag(), fma(a.imag(), b.real(), c.imag()))
		};
	}

	/** \brief Fused multiply-substract for complex&real numbers
	 * 	\return a*b+c
	*/
	template <typename number>
	inline complex<number> fma(std::complex<number> a, std::complex<number> b, number c){
		return {std::fma(a.real(),b.real(),std::fma(-b.imag(),a.imag(),c)),fms(a.real(),b.imag(),-b.real(),a.imag())};
	}

	/** \brief Fused multiply-substract for complex&real numbers
	 * 	\return a*b+c
	*/
	template <typename number>
	inline complex<number> fma(complex<number> a, number b, complex<number> c){
		return {std::fma(a.real(),b,c.real()),std::fma(a.imag(),b,c.imag())};
	}

	/** \brief Fused multiply-substract for complex&real numbers
	 * 	\return a*b+c
	*/
	template <typename number>
	inline complex<number> fma(number a,number b, complex<number> c){
		return {std::fma(a,b,c.real()),c.imag()};
	}

	/** \brief Fused multiply-substract for complex&real numbers
	 * 	\return a*b+c
	*/
	template <typename number>
	inline complex<number> fma(complex<number> a,number b, number c){
		return {std::fma(a.real(),b,c),a.imag()*b};
	}

	/** \brief Имплементация 'On the Cost of Floating-Point Computation Without Extra-Precise Arithmetic' \n
    Взято из https://github.com/KMBO19MCCO/qdrtcsKahanWepa \n
	\author Prof. W. Kahan

	Статья: https://people.eecs.berkeley.edu/~wkahan/Qdrtcs.pdf \n
	*/
	template <typename number>
    void KahanQuadratic(number a, number b, number c, vector<complex<number>> &roots){
        // x^2, x, c
        //Coefficients should be in CBA order
        try{
            b = b/static_cast<number>(-2);
            if(!std::isfinite(a)) throw std::invalid_argument("Коэффициент при x^2 равен нулю или б/м.");
            if(!std::isfinite(b)) throw std::invalid_argument("Коэффициент при x равен нулю или б/м.");
            number p = b*b;
            number q = a*c;
            //Use the hardware's FMA
            number dp = fma(b,b,-p);
            number dq = fma(a,c,-q);
            // дискриминант
            number d = (p-q) + (dp - dq);
            d = std::max(d,static_cast<number>(0));
            number S = b;
            S = std::sqrt(d)*(sign(S) + (S==0)) + S;
            number Z1 = S/a;
            number Z2 = c/S;
            if(anynotfinite(Z1, Z2)) throw std::invalid_argument("Полученные корни не определены.");

            roots = {Z1, Z2};
        }
        catch(const std::invalid_argument &err){
            std::cerr << "Error occured while working with " << a << " " << b << " " << c << " " << "\n";
            std::cerr << "Invalid argument was passed: " << err.what();
        }
        catch(const std::out_of_range &err){
            std::cerr << "Error occured while working with " << a << " " << b << " " << c << " " << "\n";
            std::cerr << "Out of range: " << err.what();
        }
        catch(...){
            std::cerr << "!!! Error occured while working with " << a << " " << b << " " << c << " " << "\n";
        }
    }

	/** \brief Случай для одного корня \n
	*/
	template <typename number>
    void simpleEquation(number a, number b, vector<complex<number>> &roots){
        // x, c
        number res = b/a;
        if(anynotfinite(a, b, res)){
            std::cerr << "Error occured while working with " << a << " " << b << " " << res << "\n";
            std::cerr << "Invalid argument was passed\n";
            return;
        }
        roots = {res};
    }

	/** \class Baydoun
	Имплементация 'Analytical formula for the roots of the general complex cubic polynomial' \n
	\author Ibrahim Baydoun 

	Статья: https://arxiv.org/abs/1512.07585 \n
	Пример использования:

		Baydoun<fp_t> Solver; // - Создаст объект Solver, который будет принимать коэффициенты указанного типа fp_t.
		int cnt = Solver(coefficients, roots); // Решит уравнение с коэффициентами, указанными в векторе coefficients,  roots - ссылка на вектор, в котором будут сохранены полученные корни.
	
	Функтор вернет int, в котором будет указано количество найденных корней.
	*/
	template <typename number>
	class Baydoun
	{
		// const number PI = static_cast<number>(_PI);
        const number PI = std::numbers::pi_v<number>;
		number onethree;
		number one27;
		number sqrt3;
		number cbrt4ftwo;

		/** \brief Аргумент комплексного числа. std::arg охватывает не все случаи
		 * 	\param inp Комплексное число.
		 * 	\param c Массив, в котором хранятся коэффициенты x.
		 * 	\param d Массив, в котором хранятся коэффициенты C.
		 * 	\return Аргумент комплексного числа.
		*/
		inline number argp(complex<number> inp){
			number x = std::real(inp);
			number y = std::imag(inp);
			if(x > 0) return std::atan2(y, x);
			else{
				return x == 0 ? (y < 0 ? -PI: PI)/static_cast<number>(2) : std::atan2(y, x) + (y < 0 ? -PI: PI);
			}
		}

		/** \brief Начало вычисления.
		 * 	\param b Массив, в котором хранятся коэффициенты x^2
		 * 	\param c Массив, в котором хранятся коэффициенты x.
		 * 	\param d Массив, в котором хранятся коэффициенты C.
		 * 	\return Количество корней.
		*/
		inline void solve(number b[6], number c[4], number d[3], vector<complex<number>> &roots){
			number d0 = d[0];
			number b0 = b[0];
			number c0 = c[0];

			number tmp1 = std::fma(static_cast<number>(2)*b0,c0,static_cast<number>(-3)*d0);
			// На основе o и r выявляется тип многочлена
			number o = std::fma(-b[1],std::fma(static_cast<number>(4)*b0,d0,-c[1]),fms(static_cast<number>(9)*d0,tmp1, static_cast<number>(4),c[2]));
			number r = std::fma(one27,fms(static_cast<number>(2), b[2], static_cast<number>(9), b0*c0),d0);
			number bthree = b0*onethree;

			// Вырожденный случай
			if(o == 0 && r == 0){
				roots = vector<complex<number>>{-bthree, -bthree, -bthree};
			}
			else{
				// Самые часто вызываемые степени
				number b0 = b[0];
				number b1 = b[1];
				number b2 = b[2];

				number c0 = c[0];
				number c1 = c[1];
				number c2 = c[2];

				number d0 = d[0];
				number d1 = d[1];
				
				number tmp = std::fma(-b0,c0,d0);

				number b0c0 = b0*c0;
				number b0c1 = b0*c1;
				number b1c1 = b1*c1;

				// Связано с комплексным многолченом в терминах его коэффициентов
				number t = fma(static_cast<number>(2)*c2, fma(static_cast<number>(33)*b0,fms(static_cast<number>(4)*b1, d0,-c0, fma(static_cast<number>(-2), d0, b0c0)), fma(static_cast<number>(4),
					fms(static_cast<number>(2),b[5],static_cast<number>(-9),d1),c2)), fma(d[2],(fma(static_cast<number>(-2),b2,static_cast<number>(9)* \
					fms(static_cast<number>(16), b0c0,static_cast<number>(3), d0))), fms(static_cast<number>(12)*b[3]*c0, fma(-static_cast<number>(7),c2,d1), b1c1*d0,fms(static_cast<number>(24),b2,static_cast<number>(-291),d0))));
				number partiond0 = fma(tmp,fms(static_cast<number>(14)*b0,c1,static_cast<number>(4)*b2,c0),fms(static_cast<number>(14)*b0c1,d0,static_cast<number>(12)*c0,d1)) + fma(b1,d1,c[3]);
				
				complex<number> sqrt1 = o >= 0 ? sqrt(o) : complex<number>(0, 1)*static_cast<number>(sqrt(fabs(o)));
				complex<number> sqrt2 = std::complex<number>(0, 1)*sqrt3;
				complex<number> sqrt2div3 = sqrt2*onethree;
				complex<number> sqrt2div9 = sqrt2div3*onethree;

				complex<number> bl = fma(sqrt1,tmp * fma(static_cast<number>(4)*b0c0,-tmp, fma(static_cast<number>(2),c2,d1)), sqrt2div9*t);
				complex<number> bl1 = pow(bl, onethree);
				complex<number> bl2 = pow(bl1, 2);
				complex<number> A1 = fma(-sqrt2div3,static_cast<number>(2)*fma(b1, fma(static_cast<number>(2),d0,tmp1),fms(static_cast<number>(15), c0*d0, static_cast<number>(13),b0c1)),static_cast<number>(2)*c0*sqrt1); 
				complex<number> A2 = fma(- sqrt1,sqrt2 *fma(static_cast<number>(2)*b0c0,fms(static_cast<number>(4),b0c0,static_cast<number>(5),d0), fma(static_cast<number>(3),d1,c2)),fma(static_cast<number>(1),
					fms(static_cast<number>(2)*b1*d0, fms(b0, d0, static_cast<number>(-58), c1), static_cast<number>(8)*b2,fms(b0c0, tmp, static_cast<number>(-5), c2)),fms(b0c0, fms(static_cast<number>(23), c2, static_cast<number>(99), d1),
					static_cast<number>(-3)*d0, fms(static_cast<number>(9), d1, static_cast<number>(7),c2))));
				complex<number>  Rbase = sqrt1 * sqrt2div9;
				complex<number> R1, R2;
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
				complex<number>  sqrt205=sqrt2*static_cast<number>(0.5);
				complex<number> M[2] = {static_cast<number>(0.5) - sqrt205, static_cast<number>(0.5) + sqrt205};

				complex<number>  arg1_1 = A1*bl1;
				complex<number>  arg1_2 = -partiond0*R1;
				complex<number>  arg2_1 = A2*bl2;
				complex<number>  arg2_2 = static_cast<number>(pow(partiond0, 2))*R2;
				
				// Вычисляем аргумент комплексного числа
				number phi1 = argp(arg1_1) - argp(arg1_2);
				number phi2 = argp(arg2_1) - argp(arg2_2); 
				complex<number> a1 = (cbrt4ftwo)*(std::cos(phi1)+std::complex<number>(0, std::sin(phi1)));
				complex<number> a2 = (cbrt4ftwo)*(std::cos(phi2)+std::complex<number>(0, std::sin(phi2)));
				roots = {
					fma(a2,R2,-fma(a1,R1,bthree)),
					fma(M[0]*a1,R1,-fma(M[1]*a2,R2,bthree)),
					fma(M[1]*a1,R1,-fma(M[0]*a2,R2,bthree))
				};
			};
		}

	public:
		Baydoun(){
			long double _onethree = 1.0L/3.0L;

			onethree = static_cast<number>(_onethree);
			one27 = static_cast<number>(pow(_onethree, 3));
			sqrt3 = static_cast<number>(std::sqrt(3L));
			cbrt4ftwo = static_cast<number>(pow(4L, _onethree)*static_cast<number>(0.5));
		}

		/** \brief Функтор для решения уравнения методом Baydoun.
		 * 	\param a Коэффициент x^3.
		 * 	\param b Коэффициент x^2.
		 * 	\param c Коэффициент x.
		 *	\param d Коэффициент C.
		 * 	\param root: Вектор, который хранит корни уравнения.
		 */
		void operator()(number a, number b, number c, number d,
				vector<complex<number>> &roots){
			// x^3, x^2, x, c
			number _b[6];
			number _c[4];
			number _d[3];
			try{
				b /= a;
				c /= a;
				d /= a;
				a = 1;
				if(anynotfinite(b, c, d)) throw std::invalid_argument("Коэффициент при x^3 равен нулю или б/м.");
				_b[0] = b; _c[0] = c; _d[0] = d;
				number bcopy = b;
				number ccopy = c;
				number dcopy = d;

				_b[1] = pow(b, 2);
				_b[2] = pow(b, 3);
				_b[3] = pow(b, 4);
				_b[4] = pow(b, 5);
				_b[5] = pow(b, 6);

				_c[1] = pow(c, 2);
				_c[2] = pow(c, 3);
				_c[3] = pow(c, 4);

				_d[1] = pow(d, 2);
				_d[2] = pow(d, 3);

				solve(_b, _c, _d, roots);
			}
			catch(const std::invalid_argument &err){
				std::cerr << "Error occured while working with " << a << " " << b << " " << c << " " << d << "\n";
				std::cerr << "Invalid argument was passed: " << err.what();
			}
			catch(const std::out_of_range &err){
				std::cerr << "Error occured while working with " << a << " " << b << " " << c << " " << d << "\n";
				std::cerr << "Out of range: " << err.what();
			}
			catch(...){
				std::cerr << "!!! Error occured while working with " << a << " " << b << " " << c << " " << d << "\n";
			}
		}

		/** \brief Функтор для решения уравнения методом Baydoun.
		 * 	\param inp Коэффициенты.
		 * 	\param root Вектор, который хранит корни уравнения.
		 * 	\param reverse Необходимо ли рассматривать коэффициенты в обратном порядке
		*/
		void operator()(vector<number> &inp, std::vector<complex<number>> &roots, bool reverse=false){
		    if(inp.size() == 2){
			reverse ?
			    simpleEquation(inp[1], inp[0], roots):
			    simpleEquation(inp[0], inp[1], roots);
		    }
		    else if(inp.size() == 3){
			reverse ?
			    KahanQuadratic(inp[2], inp[1], inp[0], roots):
			    KahanQuadratic(inp[0], inp[1], inp[2], roots);
		    }
		    else{
			reverse ?
			    operator()(inp[3], inp[2], inp[1], inp[0], roots):
			    operator()(inp[0], inp[1], inp[2], inp[3], roots);
		    }
		}
	};

	/** \class Vieta
	Имплементация Тригонометрической формулы Виеты \n
	http://poivs.tsput.ru/ru/Math/Functions/Polynomials/VietaTrigonometricFormula\n
	Пример использования:

		Vieta<fp_t> Solver; // - Создаст объект Solver, который будет принимать коэффициенты указанного типа fp_t.
		int cnt = Solver(coefficients, roots); // Решит уравнение с коэффициентами, указанными в векторе coefficients,  roots - ссылка на вектор, в котором будут сохранены полученные корни.
	
	Функтор вернет int, в котором будет указано количество найденных корней.
	*/
	template <typename number>
	class Vieta
	{
		// const number PI = static_cast<number>(_PI);
		const number PI = std::numbers::pi_v<number>;
		number sqrt3;
		number onethree;


		/** \brief Вырожденный случай.
		 * 	\param R Вычисленное значение R.
		 * 	\param b Коэффициент x^2
		 * 	\return Вектор, хранящий корни уравнения.
		*/
		inline vector<complex<number>> Degenerate(number R, number b, number inp2three){
			vector<complex<number>> roots;
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

		/** \brief Действительные корни.
		 * 	\param Q Вычисленное значение Q.
		 *	\param Q3 Вычисленное значение Q^3.
		 *	\param R Вычисленное значение R.
		 * 	\param b Коэффициент x^2
		 * 	\return Вектор, хранящий корни уравнения.
		*/
		inline vector<complex<number>> Usual(number Q, number Q3, number R, number b, number inp2three){
			vector<complex<number>> roots;
			number x1,x2,x3 = 0;


			auto _acosarg = R/sqrt(Q3);
			if(fabs(_acosarg) > 1) throw std::invalid_argument("Вызвана функция для действительных корней, но в acos |R/sqrt(Q^3)|>1!");
			number phi = acos(R/sqrt(Q3))*onethree;
			number sqrtQ = static_cast<number>(-2)*sqrt(Q);
			number half = static_cast<number>(0.5);
			number _cos = cos(phi);
			number _sin = sin(phi);

			x1 = fma(sqrtQ, _cos, -inp2three);
			x2 = fma(sqrtQ, fms(-sqrt3*half, _sin, half, _cos), -inp2three); // сos(a+b) = cos(a) * cos(b)- sin(a)*sin(b)
			x3 = fma(sqrtQ, fms(sqrt3*half,  _sin, half, _cos), -inp2three); // cos(a-b) = cos(a) * cos(b) + sin(a) * sin(b)
			roots = {x1, x2, x3};
			for (auto &r: roots) {
				if(fabs(r.imag()) < fabs(r)*std::numeric_limits<number>::epsilon()) r.imag(0);
			}
			return roots;
		}

		/** \brief Комплексные корни.
		 * 	\param Q Вычисленное значение Q.
		 *	\param Q3 Вычисленное значение Q^3.
		 *	\param R Вычисленное значение R.
		 * 	\param b Коэффициент x^2
		 * 	\return Вектор, хранящий корни уравнения.
		*/
		inline vector<complex<number>> Complex(number Q, number Q3, number R, number b, number inp2three){
			vector<complex<number>> roots;
			number x1 = 0;
			complex<number>  x2, x3 = 0;
			number _phi= 0;
			number T;
			number Tin;
			number sqrtsh;
			number sqrtabsQ3 = sqrt(fabs(Q3));
			number sqrtabsQ = sqrt(fabs(Q));

			if(sqrtabsQ3 == 0 || !std::isfinite(sqrtabsQ3)){
				return {-inp2three, -inp2three, -inp2three};
			}
			number _aharg = fabs(R)/sqrtabsQ3;
			if(Q > 0){
				if(_aharg < 1) throw std::invalid_argument("Вызвана функция для комплексных корней, но в acosh |R/sqrt(Q^3)|<1!");
				number phi = acosh(_aharg)*onethree;
				T = sqrtabsQ*cosh(phi);
				sqrtsh = sqrt3*sqrtabsQ*sinh(phi);
			}
			else{
				number phi = asinh(_aharg)*onethree;
				T = sqrtabsQ*sinh(phi);
				sqrtsh = sqrt3*sqrtabsQ*cosh(phi);
			}
			Tin = T - inp2three;
			x1 = -std::fma(static_cast<number>(2),T,inp2three);
			x2 = complex<number>(Tin,sqrtsh);
			x3 = complex<number>(Tin,-sqrtsh);
			roots = {x1, x2, x3};
			for (auto &r: roots) {
				if(fabs(r.imag()) < fabs(r)*std::numeric_limits<number>::epsilon()) r.imag(0);
			}
			return roots;
		}
	public:
		Vieta(){
			number sqrt3 = sqrt(3.0);
			number onethree = 1.0/3.0;
		}

		/** \brief Функтор для решения уравнения методом Виета.
		 * 	\param a Коэффициент x^3.
		 * 	\param b Коэффициент x^2.
		 * 	\param c Коэффициент x.
		 * 	\param d Коэффициент C.
		 * 	\param root: Вектор, который хранит корни уравнения.
		*/
		void operator()(number a, number b, number c, number d,
			vector<complex<number>> &roots){
			// x^3, x^2, x, c
			try{
				b /= a;
				c /= a;
				d /= a;
				a = 1;
				if(anynotfinite(b, c, d)) throw std::invalid_argument("Коэффициент при x^3 равен нулю или б/м.");

				number b_onethree = b*onethree;
				auto Q = fms(b_onethree,b_onethree,c,onethree);
				if(Q == 0 || !std::isfinite(Q)){
					number rs = -b_onethree;
					roots = {rs, rs, rs};
				}
				else{
					number R = std::fma(static_cast<number>(0.5), std::fma(-c, b_onethree, d), pow(b_onethree, 3));
					auto Q3 = pow(Q, 3);
					auto S = std::fma(-R,R,Q3);
					if(S==0 || !std::isfinite(S)){
						roots = Degenerate(R, b, b_onethree);
					}
					else if(S > 0){
						roots = Usual(Q, Q3, R, b, b_onethree);
					}
					else{
						roots = Complex(Q, Q3, R, b, b_onethree);
					}
				}
			}
			catch(const std::invalid_argument &err){
				std::cerr << "Error occured while working with " << a << " " << b << " " << c << " " << d << "\n";
				std::cerr << "Invalid argument was passed: " << err.what();
			}
			catch(const std::out_of_range &err){
				std::cerr << "Error occured while working with " << a << " " << b << " " << c << " " << d << "\n";
				std::cerr << "Out of range: " << err.what();
			}
			catch(...){
				std::cerr << "!!! Error occured while working with " << a << " " << b << " " << c << " " << d << "\n";
			}
		}

		/** \brief Функтор для решения уравнения методом Виета.
		 * 	\param inp Коэффициенты.
		 * 	\param root Вектор, который хранит корни уравнения.
		 * 	\param reverse Необходимо ли рассматривать коэффициенты в обратном порядке
		*/
		void operator()(vector<number> &inp, std::vector<complex<number>> &roots, bool reverse=false){
		    if(inp.size() == 2){
			reverse ?
			    simpleEquation(inp[1], inp[0], roots):
			    simpleEquation(inp[0], inp[1], roots);
		    }
		    else if(inp.size() == 3){
			reverse ?
			    KahanQuadratic(inp[2], inp[1], inp[0], roots):
			    KahanQuadratic(inp[0], inp[1], inp[2], roots);
		    }
		    else{
			reverse ?
			    operator()(inp[3], inp[2], inp[1], inp[0], roots):
			    operator()(inp[0], inp[1], inp[2], inp[3], roots);
		    }
	}
};
}
#endif
