#include <stdio.h>

#include <stdlib.h>
#include <math.h>


double sgn(double x) { // Signumfunction except sign(0) = 1
	return x >= 0 ? 1.0 : -1.0;
}

double function_f(double x) {
	return exp(-x*x);
}

double function_g(double x) {
	return 1.0/((x+1)*sqrt(x));
}

double function_h(double x) {
	return sin(sqrt(x)) * exp(sqrt(x)) / sqrt(x);
}

double eval_improper_integral(double func (double), double a, double b, double precision) {
	double sum = 0.0;
	if(isinf(a) && isinf(b)) { 					// Assuming convergence and diagonal and vertical symmetry
		a = 0.0; 								// Variable a is reused
		double delta = precision, py = func(a), ny = func(-a); 	// positive y, negative y
		do {
			if(isfinite(py) && py != 0) { // symmetry assumption
				delta = fabs(precision/py);
				sum += sgn(py) * precision;
				sum += sgn(ny) * precision;
			}
			a+=delta;
			py = func(a);
			ny = func(-a);
		} while(1.0/delta > precision);
	} else if(isinf(b)) {
		double delta = precision, y = func(a); 	// positive y, negative y
		do {
			if(isfinite(y) && y != 0) {
				delta = fabs(precision/y);
				sum += sgn(y) * precision;
			}
			a+=delta;
			y = func(a);
		} while(1.0/delta > precision);
	}
	return sum;
}

double eval_integral(double func (double), double a, double b, double precision) {
	double sum = 0.0, k = 1.0;
	if(a > b) {
		int swap = a;
		a=b;
		b=swap;
		k = -k;
	}
	if(isfinite(a) && isfinite(b)) {
		double delta = precision, y;
		while(a <= b) {
			y = func(a);
			if(isfinite(y) && y != 0) {
				delta = fabs(precision/y); // Adjust delta accordingly
				sum += sgn(y)*precision; // y*delta = y*p1/y = p1
			}
			a+=delta;
		}
	} else {
		return eval_improper_integral(func, a, b, precision);
	}
	return k * sum;
}

int main() {
	printf("a=0, b=2, sin(sqrt(x))e^sqrt(x))/sqrt(x) ≈ %.6f\n", eval_integral(&function_h, 0, 2, 0.000001));
	printf("a=-∞, b=∞, ∫e^{-x^2} ≈ %.6f\n", eval_integral(&function_f, -INFINITY, INFINITY, 0.0000001));
	printf("a=0, b=∞, ∫1/{(x+1)sqrt(x)} ≈ %.6f\n", eval_integral(&function_g, 0, INFINITY, 0.0000001));
}
