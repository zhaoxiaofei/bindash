#include "../bhmath.hpp"

#include <stdio.h>

int test_pois() {
	double res1 = bhmath_ppois_tail(0, 5);
	double res2 = bhmath_ppois_tail(1, 5);
	double res3 = bhmath_ppois_tail(5, 5);
	double res4 = bhmath_ppois_tail(10, 5);
	double res5 = bhmath_ppois_tail(20, 5);
	double res6 = bhmath_ppois_tail(40, 5);
	printf("ppois_tail(0, 5) = %.9e\n", res1);
	printf("ppois_tail(1, 5) = %.9e\n", res2);
	printf("ppois_tail(5, 5) = %.9e\n", res3);
	printf("ppois_tail(10, 5) = %.9e\n", res4);
	printf("ppois_tail(20, 5) = %.9e\n", res5);
	printf("ppois_tail(40, 5) = %.9e\n", res6);
	printf("ppois_tail(40, 3) = %.9e\n", bhmath_ppois_tail(40, 3));
}

int test_norm(int tablesize = 300) {
	for (int i = 0; i < tablesize; i+=20) {
		printf("z=%lf,1-cdf(z)=%.9e\n", (double)i / 20, NORMTABLE20TAIL[i]); // i * 10 / tablesize
	}
}

void test_binom_tail_print(unsigned int x,  unsigned int n, double p) {
	double res = bhmath_pbinom_tail(x, n, p);
	printf("pbinom_tail(%u, %u, %f)==%.9e\n", x, n, p, res);
}

int test_binom(int tablesize = 300) {
	
	test_binom_tail_print(  0,  100, 0.1);
	test_binom_tail_print(  1,  100, 0.1);
	test_binom_tail_print( 10,  100, 0.1);
	test_binom_tail_print( 90,  100, 0.1);
	test_binom_tail_print( 99,  100, 0.1);
	test_binom_tail_print(100,  100, 0.1);
	
	test_binom_tail_print(  0, 1000, 0.0049);
	test_binom_tail_print(  0, 1000, 0.0051);
	test_binom_tail_print(  0, 1001, 0.0049);
	test_binom_tail_print(  0, 1001, 0.0051);
	test_binom_tail_print(  1, 1000, 0.0049);
	test_binom_tail_print(  1, 1000, 0.0051);
	test_binom_tail_print(  1, 1001, 0.0049);
	test_binom_tail_print(  1, 1001, 0.0051);

	test_binom_tail_print(  5, 1000, 0.0049);
	test_binom_tail_print(  5, 1000, 0.0051);
	test_binom_tail_print( 10, 1000, 0.0049);
	test_binom_tail_print( 10, 1000, 0.0051);
	test_binom_tail_print( 20, 1000, 0.0049);
	test_binom_tail_print( 20, 1000, 0.0051);
}

int main(int argc, char **argv) {
	test_pois();
	test_norm();
	test_binom();
}

