#include "SMO_LIMITED.h"
#include <stdio.h>
#include <stdlib.h>

void SMO_LIMITED :: plot_len_time_work(int n, FILE* p, FILE* f0, FILE* f1, FILE* f2, FILE* f3){
	double p0 = 0.0;
	double pnm = 1.0;
	double pn = 1.0;
	double sum = 0.0;
	int m;

	for(int i = n; i >= 1; i--){
		m = n - i;
		for(int k = 0; k <= i; k++){
			sum += pow(a, k) / fact(k);
		}
		sum += pow(a, i + 1) / (i * fact(i)) * (1.0 - pow(a / i, m)) / (1.0 - a / i);

		p0 = 1.0 / sum;
		//pn = p0 * pow(a, i) / fact(i);
		pnm = p0 * pow(a, i + m) / (pow(i, m) * fact(i));

		sum = 0.0;

		double K = a * (1.0 - pnm) / i;
		double Len = pow(a, i + 1) * p0 / (fact(i - 1) * pow(i - a, 2.0));
		double L = ((pow(a, i + 1)*p0) / (i*fact(i))*((1.0 - pow(a / i, m)*(m + 1.0 - m*a / i)) / pow(1.0 - a / i, 2)));
		double T = L * tc;

		/*printf("LTW %d %.10lf\n", i, pnm);
		printf("LTW %d %.10lf\n", i, K);
		printf("LTW %d %.10lf\n", i, Len);
		printf("LTW %d %.10lf\n", i, L);
		printf("LTW %d %.10lf\n", i, T);

		getchar();*/

		fprintf(p, "%d %.10lf\n", i, pnm);
		fprintf(f0, "%d %.10lf\n", i, K);
		if(i >= a + 1){
		fprintf(f1, "%d %.10lf\n", i, Len);
		fprintf(f2, "%d %.10lf\n", i, L/i);
		fprintf(f3, "%d %.10lf\n", i, T);
		}
	}	
}

void SMO_LIMITED :: find_serv_queue(){
	int n, m;
	double p0 = 0.0;
	double pnm = 1.0;
	double sum = 0.0;
	double eps = 1.0e-3;
	int max = 14;
	double L[max][max] = {0.0};
	double K[max][max] = {0.0};

	for(n = a + 1; n < max; n++){
		for(m = 1; m < max; m++){
			for(int k = 0; k <= n; k++){
				sum += pow(a, k) / fact(k);
			}
			sum += pow(a, n + 1) / (n * fact(n)) * (1.0 - pow(a / n, m)) / (1.0 - a / n);
			p0 = 1.0 / sum;
			pnm = p0 * pow(a, n + m) / (pow(n, m) * fact(n));
			sum = 0.0;

			K[n][m] = a * (1.0 - pnm);
			L[n][m] = ((pow(a, n + 1)*p0) / (n*fact(n))*((1.0 - pow(a / n, m)*(m + 1.0 - m*a / n)) / pow(1.0 - a / n, 2)));
		}
	}

	double prev = 0.0;
	double cur = 0.0;
	int ind = a + 1, jnd = 1;

	for(int i = a + 1; i < max; i++){
		cur = fabs(K[i][1] / i  - L[i][1] / 1);
		prev = cur;
		for(int j = 1; j < max; j++){
			cur = fabs(K[i][j] / i - L[i][j] / j);
			printf("K %lf L %lf R %lf C %lf P %lfn = %d m = %d", K[i][j] / i, L[i][j] / j, fabs(K[i][j] / i - L[i][j] / j), cur, prev, i, j);
			getchar();
			if(cur > prev){
				ind = i;
				jnd = j - 1;
				printf("3 is === %d %d %lf\n", ind, jnd, prev);
				printf("%lf %lf\n", K[ind][jnd] / ind, L[ind][jnd] / jnd);
				return;
			}
			prev = cur;
		}
	}
}
