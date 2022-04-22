#include <stdio.h>
#include <math.h>

double function_f(double x) {
	return 1.0f/((x+1)*sqrt(x));
	//return exp(-x*x);
}

double solve_integral(double func (double), double a, double b) {
	const double p1 = 0.000001f, p2 = 0.000002f;
	double sum = 0.0f, delta = p1, temp;
	while(a <= b) {
		temp = func(a);
		if(isfinite(temp)) {
			if(temp*delta > p2 || temp*delta < p1) delta = p1/temp;
			sum += func(a)*delta;
		}
		a+=delta;
	}
	return sum;
}

int main() {
	double a = 0.0f;
	double b = 100000.0f;

	printf("Evaluation: %.15f", solve_integral(&function_f, a, b));
}
