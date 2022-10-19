#ifndef BAYDOUN_H
#define BAYDOUN_H
#include <cmath>
#include <vector>
#include <complex>
#include <iostream>

using namespace std::complex_literals;



template <typename number>
class Baydoun
{
	number onethree;
	number twothree;
	number one27;

	number sqrt3;
	number cbrt4;

	number bthree;

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

	std::vector<std::complex<number>> part2(number *a, number *b, number *c, number *d,
		number o, number r){
		number c0d0 = c[0]*d[0];
		number b0c0 = b[0]*c[0];
		number b0c1 = b[0]*c[1];
		number b1c1 = b[1]*c[1]*4;
		number t = c[2]*(16*b[5] + 72*d[1] + 264*b[2]*d[0]  + 66*b[1]*c[1] - \
		    132*b[0]*c0d0 + 2*c[2]) + b[3]*c[0]*(12*d[1] - 84*c[2]) - b[1] * \
		    c[1]*d[0]*(24*b[2]+291*d[0]) + d[2]*(144*b0c0 - 27*d[0] - 2*b[2]);
		number d0 = 4*(b[3]*c[1] - b[2]*c0d0) - 12*c[0]*d[1] -14*b[1]*c[2] + \
			28*b0c1*d[0] + b[1]*d[1] + c[3];

		std::cout << "_P2 c0, b0c0, b0c1, b1c1, t, d0 " << c0d0 << " " << b0c0 << 
			" " << b0c1 << " " << b1c1 << " " << t << " " << d0 << "\n";

		std::complex<number> sqrt1;
		if (o > 0)
		    sqrt1 = pow(std::complex<number>(o, 0), 0.5);
		else
		    sqrt1 = std::complex<number>(0, 1)*pow(std::complex<number>(abs(o), 0), 0.5);
		std::cout << "_P2 sqrt1 " << sqrt1 << "\n";
		std::complex<number> sqrt2 = std::complex<number>(0, 1)*sqrt3;
		std::complex<number> sqrt3 = cbrt4;

		auto sqrt2div3 = sqrt2*onethree;
		auto sqrt2div9 = sqrt2*onethree*onethree;
		auto sqrt3ftwo = sqrt3*0.5;

		auto bl = (d[0]-b0c0) * sqrt1 * (b1c1 - 4*b[0]*c0d0 + 2*c[2] + d[1]) +sqrt2div9*t;
		auto bl1 = pow(bl, onethree);
		auto bl2 = pow(bl1, 2);
		auto A1 = (-sqrt2div3)*(8*b[2]*c[0] - 4*d[0]*b[1] - 26*b0c1 + 30*c0d0) + 2*c[0]*sqrt1;
		auto A2 = 8*(b[4]*c[1] - b[3]*c0d0) - 40*b[2]*c[2] + 2*b[2]*d[1] +\
		    29*b1c1*d[0] + 23*b[0]*c[3] -21*c[2]*d[0] +\
		    27*d[2] - 99*b0c0*d[1] -sqrt1*sqrt2 * (2*b1c1 - 10*b[0]*c0d0 + c[2] + 3*d[1]);
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


		auto sqrt205=sqrt2*0.5;
		std::complex<number> M[2] = {0.5 - sqrt205, 0.5 + sqrt205};
		std::complex<number> M2[2] = {-0.5 - sqrt205, -0.5 + sqrt205};

		auto arg1_1 = A1*bl1;
		auto arg1_2 = -d0*R1;
		auto arg2_1 = A2*bl2;
		auto arg2_2 = pow(d0, 2)*R2;

		// Вычисляем аргумент комплексного числа
		auto phi1 = std::arg(arg1_1) - std::arg(arg1_2);
		auto phi2 = std::arg(arg2_1) - std::arg(arg2_2);

		auto a1 = (sqrt3ftwo)*(std::cos(phi1)+1i*std::sin(phi1));
		auto a2 = (sqrt3ftwo)*(std::cos(phi2)+1i*std::sin(phi2));

		auto a1R1 = a1*R1;
		auto a2R2 = a2*R2;
		auto x1 = -a1R1 + a2R2 - bthree;
		auto x2 = M[0]*a1R1 + M2[0]*a2R2 - bthree;
		auto x3 = M[1]*a1R1 + M2[1]*a2R2 - bthree;

		return std::vector<std::complex<number>>{x1, x2, x3};
	}

	std::vector<std::complex<number>> solve(number *a, number *b, number *c, number *d){
		bthree = b[0]*onethree;
		number o = -4*(b[2]*d[0] + c[2]) + b[1]*c[1] + 18*b[0]*c[0]*d[0] - 27*d[1];
		number r = 2*b[2]*one27 -9*b[0]*c[0]*one27 + d[0];

		std::cout << "_SOLVE o r, one27" << o << " " << r << " " << one27 << "\n";

		if(o == 0 && r == 0){
		    return std::vector<std::complex<number>>{-bthree, -bthree, -bthree};
		}
		else{
		    return part2(a, b, c, d, o, r);
		};
	}

public:

	Baydoun(){
		onethree = 1.0/3.0;
		twothree = 2*onethree;
		one27 = onethree*onethree*onethree;

		sqrt3 = std::sqrt(3);
		cbrt4 = std::cbrt(4);
	}

	int operator()(number a, number b, number c, number d,
			std::vector<std::complex<number>> &roots){
			// x^3, x^2, x, c
		std::cout << "onethree, twothree, one27, sqrt3, cbrt4 " <<
			onethree << " " << twothree << " " << one27 << " " << sqrt3 << " " << cbrt4 << "\n";
		std::cout << a << " " << b << " " << c << " " << d << "\n";
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
		number *_a = new number[1];
		_a[0] = a; _b[0] = b; _c [0]= c; _d[0] = d;
		prepare(_b, _c, _d);

		std::cout << "\nDone preparing\n";
		auto result = solve(_a, _b, _c, _d);
		std::cout << "Done solving\n";
		roots = result;
		return 0;
	}
};

#endif
