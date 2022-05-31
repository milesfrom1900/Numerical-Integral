//By Sebastian Miles 2022

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
	//return 1.0/((x+1)*sqrt(x));
	//return 2*exp(x)/((exp(x)-1)*(exp(x)-1)+1);
	//return 2*2*pow(x,1)/(pow(x,4)+1);
	return sqrt(1-x*x);
}

double function_h(double x) {
	return sin(sqrt(x)) * exp(sqrt(x)) / sqrt(x);
}

double function_a(double x) {
	return x;
}

double eval_improper_integral(double func (double), double a, double b, double precision) {
	double sum = 0.0;
	if(isinf(a) && isinf(b)) { 					// Assuming convergence and diagonal and vertical symmetry
		a = 0.0; 								// Variable a is reused
		double delta = 0.00000000000001, py = func(a), ny = func(-a), lpy = py, lny = ny; 	// positive y, negative y
		do {
			if(isfinite(py) && py != 0) { // symmetry assumption
				sum += sgn(py) * (precision - fabs(lpy-py)*0.5*delta);
				sum += sgn(ny) * (precision - fabs(lny-ny)*0.5*delta);
				//sum += sgn(py) * (precision);
				//sum += sgn(ny) * (precision);
				delta = fabs(precision/py);
			}
			a+=delta;
			lpy = py;
			py = func(a);
			lny = ny;
			ny = func(-a);
			if(!isfinite(lpy)) lpy = py;
			if(!isfinite(lny)) lny = ny;
		} while(a < 1e10);
	} else if(isinf(b)) {
		double delta = 0.00000000000001, y = func(a), ly = y;
		do {
			if(isfinite(y) && y != 0) {
				sum += sgn(y) * (precision);
				delta = fabs(precision/y);
			}
			a+=delta;
			ly = y;
			y = func(a);
			if(!isfinite(ly)) ly = y;
		} while(a < 1e10);
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
		double delta = precision, y, ly = func(a);
		while(a <= b) {
			y = func(a);
			if(isfinite(y) && y != 0) {
				sum += (y+sgn(y)*fabs(ly-y)*0.5)*delta; // y*delta = y*p1/y = p1
			}
			a+=delta;
			ly=y;
		}
	} else {
		return eval_improper_integral(func, a, b, precision);
	}
	return k * sum;
}

int main() {
	//printf("a=-∞, b=∞, ∫e^{-x^2}dx ≈ %.6f\n", eval_integral(&function_f, -INFINITY, INFINITY, 0.000002));
	printf("a=0, b=∞, ∫1/{(x+1)sqrt(x)}dx ≈ %.7f\n", 2*eval_integral(&function_g, -1, 1, 0.0000001));
	
	//printf("a=0, b=1, ∫xdx ≈ %.6f\n", eval_integral(&function_a, 0, 1, 0.00001));
}
