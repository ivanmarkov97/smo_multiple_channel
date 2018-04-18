#include <stdio.h>
#include <stdlib.h>
#include "SMO.h"

int SMO :: num_oper(double delta, FILE* f){
	double pn = 1.0;
	double p0 = 0.0;
	int n = -1;
	double sum = 0.0;
	while(pn >= delta){
		n++;
		sum += pow(a, n) / fact(n);
		//printf("sum == %lf\n", sum);
		p0 = 1 / sum; 
		pn = p0 * pow(a, n) / fact(n);
		//printf("a == %lf p0 == %lf pn == %lf, n == %d\n", pow(a, n), p0, pn, n);
		fprintf(f, "%d %.10lf\n", n, pn);
		//getchar();
	}
	return n == -1? 0 : n;
}

void SMO :: plot_load(int n, FILE* k){
	double sum = 0.0;
	double p0 = 0.0;
	double pn = 1.0;
	for(int i = 1; i <= n; i++){
		for(int j = 0; j <= i; j++){
			sum += pow(a, j) / fact(j);
		}
		p0 = 1.0 / sum; 
		pn = p0 * pow(a, i) / fact(i);
		sum = 0.0;
		fprintf(k, "%d %.10lf\n", i, a * (1.0 - pn) / i);
	}
}

/*void SMO :: plot_len_time_work(int n, FILE* p, FILE* f0, FILE* f1, FILE* f2, FILE* f3){
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

		getchar();

		fprintf(p, "%d %.10lf\n", i, pnm);
		fprintf(f0, "%d %.10lf\n", i, K);
		if(i >= a + 1){
		fprintf(f1, "%d %.10lf\n", i, Len);
		fprintf(f2, "%d %.10lf\n", i, L/i);
		fprintf(f3, "%d %.10lf\n", i, T);
		}
	}	
}*/

void SMO :: task_2(FILE* f1, FILE* f2){
	//396 11 45
	tc = 279.0;
	ts = 9.0;
	a = ts / tc;
	int N = 38;
	int m;
	double p0 = 0.0;
	double sum = 0.0;

	double PI[N + 1][N + 1] = {0.0};
	double N_sr = 0.0;
	double K_load = 0.0;

	for(int m = 1; m <= N; m++){
		for(int i = 0; i <= m; i++){
			sum += pow(a, i) * fact_rev(N, i) / fact(i);
			PI[m][i] = pow(a, i) * fact_rev(N, i) / fact(i);
		}
		for(int j = m + 1; j <= N; j++){
			sum += pow(a, j) * fact_rev(N, j) / (fact(m) * pow(m, j - m));
			PI[m][j] = pow(a, j) * fact_rev(N, j) / (fact(m) * pow(m, j - m));
		}
		p0 = 1.0 / sum;
		sum = 0.0;
		//printf("task 2 : m == %d  p0 == %.10lf\n", m, p0);

		for(int k = 0; k <= N; k++){
			N_sr += (N - k) * PI[m][k] * p0;
		}
		//printf("N_sr == %.20lf\n", N - N_sr);
		fprintf(f1, "%d %.20lf\n", m, N - N_sr);

		//printf("pn == %.20lf\n", PI[m][N]);

		for(int i = 0; i <= m; i++){
			K_load += i * p0 * PI[m][i];
		}

		for(int j = m + 1; j <= N; j++){
			K_load += m * p0 * PI[m][j];
		}

		//printf("K_load == %.20lf\n", K_load);
		fprintf(f2, "%d %.20lf\n", m, K_load / m);
		K_load = 0.0;
		N_sr = 0.0;
		//getchar();
	}
}
