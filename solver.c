#include <stdio.h>
#include <math.h>

double function_f(double x) {
	return 1.0f/((x+1)*sqrt(x));
	//return exp(-x*x);
}

double solve_integral(double func (double), double a, double b) {
	double sum = 0.0f, delta = 0.0001f, temp;
	while(a <= b) {
		temp = func(a);
		if(isfinite(temp))
			sum += temp*delta;
		a+=delta;
	}
	return sum;
}

int main() {
	double a = 0.0f;
	double b = 1000.0f;

	printf("Evaluation: %.15f", solve_integral(&function_f, a, b));
}
