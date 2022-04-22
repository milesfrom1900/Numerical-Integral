#include <stdio.h>

// abs()
#include <stdlib.h>
#include <math.h>


double sgn(double x) { // Signumfunction except sign(0) = 1
	return x >= 0 ? 1.0 : -1.0;
}

double function_f(double x) {
	//return 1.0/((x+1)*sqrt(x));
	return exp(-x*x);
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
		} while(1.0/delta > precision || py==0); // symmetry assumption
	} else if(isinf(b)) {
		double delta = precision, y = func(a); 	// positive y, negative y
		
		do {

			if(isfinite(y) && y != 0) { // symmetry assumption
				delta = fabs(precision/y);
				sum += sgn(y) * precision;
			}
			a+=delta;
			y = func(a);
			//printf("%f\n", fabs(delta*y));
			//printf("%i\n", fabs(delta*y) >= 0.01);
		} while(1.0/delta > precision || y==0); // symmetry assumption
	}
	return sum;
}

double eval_integral(double func (double), double a, double b, double precision) {
	// Precision
	const double p1 = precision;

	double sum = 0.0, k = 1.0;
	
	if(a > b) {
		int swap = a;
		a=b;
		b=swap;
		k = -k;
	}
	if(isfinite(a) && isfinite(b)) {
		double delta = p1, y;
		while(a <= b) {
			y = func(a);
			if(isfinite(y) && y != 0) {
				delta = fabs(p1/y); // Adjust delta accordingly
				sum += sgn(y)*p1; // y*delta = y*p1/y = p1
			}
			a+=delta;
		}
	} else {
		return eval_improper_integral(func, a, b, precision);
	}
	return k * sum;
}

int main() {
	double a = 0.0;
	double b = INFINITY;

	printf("Evaluation: %.15f", eval_integral(&function_f, a, b, 0.00001));
}
