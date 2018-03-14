#include "SMO_LIMITED.h"
#include <stdio.h>
#include <stdlib.h>

void SMO_LIMITED :: find_serv_queue(){
	int n, m;
	double p0 = 0.0;
	double pnm = 1.0;
	double sum = 0.0;
	double eps = 1.0e-3;
	int max = 14;

	double L[max][max] = {0.0};
	double K[max][max] = {0.0};

	for(n = 5 + 1; n < max; n++){
		for(m = 0; m < max; m++){
			for(int k = 0; k <= n; k++){
				sum += pow(a, k) / fact(k);
			}
			sum += pow(a, n + 1) / (n * fact(n)) * (1.0 - pow(a / n, m)) / (1.0 - a / n);
			p0 = 1.0 / sum;
			pnm = p0 * pow(a, n + m) / (pow(n, m) * fact(n));
			sum = 0.0;

			K[n][m] = a * (1.0 - pnm); //pQ
			L[n][m] = p0 * (pow(a, n + 1) / (n * fact(n))) * ( 1.0 - pow(a / n, m) * (m + 1.0 - m / n * a)) / pow(1.0 - a / n, 2.0);
			L[n][m] = ((pow(a, n + 1)*p0) / (n*fact(n))*((1.0 - pow(a / n, m)*(m + 1.0 - m*a / n)) / pow(1.0 - a / n, 2)));
		}
	}

	double del = fabs(K[7][0] - L[7][0]);
	int ind = 5 + 1, jnd = 0;

	for(int i = 6; i < n; i++){
		for(int j = 0; j < m; j++){

			//printf("%lf %lf %d %d\n", K[i][j], L[i][j], i, j);
			//getchar();

			if(fabs(K[i][j] - L[i][j]) < del /*eps*/){
				del = fabs(K[i][j] - L[i][j]);
				ind = i;
				jnd = j;
				//printf("K = %lf L = %lf n = %d  m = %d\n", K[i][j], L[i][j], i, j);
			}
		}
		//putchar('\n');
	}
	printf("3 is === %d %d %lf\n", ind, jnd, del);
	printf("%lf %lf\n", K[ind][jnd], L[ind][jnd]);
}
