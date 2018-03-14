#include "SMO_UNLIMITED.h"
#include <stdio.h>
#include <stdlib.h>

void SMO_UNLIMITED :: plot_len_time(int n, FILE* f0, FILE* f1, FILE* f2){
	double p0 = 0.0;
	double pn = 1.0;
	double sum = 0.0;
	double L = 0.0;
	double T = 0.0;

	for(int i = 7; i <= n; i++){
		for(int k = 0; k <= i; k++){
			sum += pow(a, k) / fact(k);
		}
		sum += pow(a, i + 1) / (fact(i) * (i - a));
		p0 = 1.0 / sum;
		sum = 0.0;

		fprintf(f0, "%d %lf\n", i, a / i);

		printf(" p0 ===== %lf %d\n", p0, i);

		L = p0 * pow(a, i + 1) / (i * fact(i) * pow(1.0 - a / i, 2));
		fprintf(f1, "%d %lf\n", i, L);
		printf("%d %lf\n", i, L);

		T = (L) * tc; 

		fprintf(f2, "%d %lf\n", i, T/i);
		printf("%d %lf\n", i, T);

		//getchar();
	}
}

double SMO_UNLIMITED :: calc_p0_5(int n, int m){
	double p0 = 0.0;
	double sum = 0.0;
	double sum_b = 0.0;
	double znam = 1.0;

	for(int i = 0; i <= n; i++){
		sum += pow(a, i) / fact(i);
	}
	for(int j = 1; j <= m; j++){
		znam *= n + j * b;
		sum_b += pow(a, j) / znam;
	}
	sum += (pow(a, n) / fact(n)) * sum_b;
	//printf(" sum_b == %.10lf\n", sum_b);
	p0 = 1.0 / sum;
	return p0;
}

void SMO_UNLIMITED :: task_exit(int n, FILE* k, FILE* l, FILE* t){
	double eps = 1.0e-15;
	double p0 = 0.0;
	double pnm = 1.0;
	double K = 0.0;
	double L = 0.0;
	double T = 0.0;
	int m = 1;

	//printf("start m = 1 %.10lf\n", calc_p0_5(tc, ts, tw, 12, 15));
	//printf("start m = 2 %.10lf\n", calc_p0_5(tc, ts, tw, n, 2));

	for(int i = 7; i <= n; i++){
		while(fabs(calc_p0_5(i, m) - calc_p0_5(i, m + 1) > eps)){
			m++;	
			p0 = calc_p0_5(i, m);
		}
		//printf("%d %d\n", i, m);
		double znam = 1.0;
		for(int j = 1; j <= m; j++){
			znam *= (i + j * b);
		}
		pnm = pow(a, i + m) / (fact(i) * znam);
		//pnm = pow(a, i + m) / (pow(i, m) * fact(i)) * p0;
		double prefix = pow(a, i) / fact(i);

		for(int j = 1; j <= m; j++){
			double Pnj = 1.0;
			for(int k = 1; k <= j; k++){
				Pnj *= a / (n + k * b);
			}
			L += j * Pnj;
		}

		L *= prefix;
		//printf("pnm == %.10lf\n", pnm);

		K = a * (1.0 - pnm) / i ;/// i;
		//L = ((pow(a, i + 1)*p0) / (i*fact(i))*((1.0 - pow(a / i, m)*(m + 1.0 - m*a / i)) / pow(1.0 - a / i, 2)));
		T = (L) * tc; 

		fprintf(k, "%d %lf\n",i, K);
		fprintf(l, "%d %lf\n",i, L);
		fprintf(t, "%d %lf\n",i, T);

		znam = 1.0;
		L = 0.0;

		//getchar();
		m = 1;
	}
}
