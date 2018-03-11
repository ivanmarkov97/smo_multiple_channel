#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#define INVALID_ARG -1.0
#define ERROR_OPEN_FILE -2.0

double fact(int n){
	double ret = 1.0;
	for(int i = 1; i <= n; i++){
		ret *= i;
	}
	return ret;
}

double fact_rev(int n, int m){
	double ret = 1.0;
	ret = fact(n) / fact(n - m);
	return ret;
}

int operators_num(double tc, double ts, double delta, FILE* f){
	double a = ts / tc;
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

void find_serv_queue(double tc, double ts){
	int n, m;
	double p0 = 0.0;
	double pnm = 1.0;
	double sum = 0.0;
	double eps = 1.0e-3;
	int max = 20;

	double L[max][max] = {0.0};
	double K[max][max] = {0.0};

	double a = ts / tc;

	for(n = 1; n < max; n++){
		for(m = 0; m < max; m++){
			for(int k = 0; k <= n; k++){
				sum += pow(a, k) / fact(k);
			}
			sum += pow(a, n + 1) / (n * fact(n)) * (1.0 - pow(a / n, m)) / (1.0 - a / n);
			p0 = 1.0 / sum;
			pnm = p0 * pow(a, n + m) / (pow(n, m) * fact(n));
			sum = 0.0;
			if(n == 5 && m == 5){
				printf("%lf\n", p0);
				//getchar();
			}
			printf("p0 == %lf pn+m == %lf n == %d m == %d\n", p0, pnm, n, m);
			//getchar();
			K[n][m] = a * (1.0 - pnm); //pQ
			L[n][m] = p0 * (pow(a, n + 1) / (n * fact(n))) * ( 1.0 - pow(a / n, m) * (m + 1.0 - m / n * a)) / pow(1.0 - a / n, 2.0);
			L[n][m] = ((pow(a, n + 1)*p0) / (n*fact(n))*((1.0 - pow(a / n, m)*(m + 1.0 - m*a / n)) / pow(1.0 - a / n, 2)));
		}
	}

	for(int i = 1; i < n; i++){
		for(int j = 0; j < m; j++){
			if(fabs(K[i][j] - L[i][j]) < eps){
				printf("K = %lf L = %lf n = %d  m = %d\n", K[i][j], L[i][j], i, j);
			}
		}
		putchar('\n');
	}

	printf("%lf %lf\n", K[3][4], L[3][4]);
}

void plot_load(double tc, double ts, int n, FILE* k){
	double a = ts / tc;
	double sum = 0.0;
	double p0 = 0.0;
	double pn = 1.0;
	for(int i = 1; i <= n; i++){
		sum += pow(a, i) / fact(i);
		p0 = 1 / sum; 
		pn = p0 * pow(a, i) / fact(i);
		fprintf(k, "%d %.10lf\n", i, a * (1.0 - pn) / n);
	}
}

void plot_reverse(double tc, double ts, int n, FILE* r, FILE* k){
	double a = ts / tc;
	double sum = 0.0;
	double p0 = 0.0;
	double pn = 1.0;
	for(int j = n; j >= 0; j--){
		for(int i = 0; i <= j; i++){
			sum += pow(a, i) / fact_rev(j, i);
		}
		p0 = 1.0 / sum;
		pn = p0 * pow(a, j) / fact(j);
		sum = 0.0;
		//printf("fact_rev == %lf sum == %lf p0 == %lf pn == %lf, n == %d\n", fact_rev(n, j), sum, p0, pn, j);
		fprintf(r, "%d %.10lf\n", j, pn); //(a* (1 - pn) / n)
		fprintf(k, "%d %.10lf\n", j, a * (1.0 - pn) / n);
		//getchar();
	}
}

void plot_len_time(double tc, double ts, int n, FILE* f1, FILE* f2){
	double p0 = 0.0;
	double pn = 1.0;
	double sum = 0.0;
	double a = ts / tc;

	double L = 0.0;
	double T = 0.0;

	for(int i = 6; i <= n; i++){
		for(int k = 0; k <= i; k++){
			sum += pow(a, k) / fact(k);
		}
		sum += pow(a, i + 1) / (fact(i) * (i - a));
		p0 = 1.0 / sum;
		sum = 0.0;

		printf(" p0 ===== %lf %d\n", p0, i);

		L = p0 * pow(a, i + 1) / (i * fact(i) * pow(1.0 - a / i, 2));
		fprintf(f1, "%d %lf\n", i, L);
		printf("%d %lf\n", i, L);

		T = (L) * tc; 

		fprintf(f2, "%d %lf\n", i, T);
		printf("%d %lf\n", i, T);

		getchar();
	}
}

double calc_p0_5(double tc, double ts, double tw, int n, int m){
	double a = ts / tc;
	double b = ts / tw;
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

void task_exit(double tc, double ts, double tw, int n, FILE* k, FILE* l, FILE* t){
	double a = ts / tc;
	double eps = 1.0e-15;

	double p0 = 0.0;
	double pnm = 1.0;
	
	double K = 0.0;
	double L = 0.0;
	double T = 0.0;

	int m = 1;

	printf("start m = 1 %.10lf\n", calc_p0_5(tc, ts, tw, 12, 15));
	printf("start m = 2 %.10lf\n", calc_p0_5(tc, ts, tw, n, 2));

	for(int i = 6; i <= n; i++){
		while(fabs(calc_p0_5(tc, ts, tw, i, m) - calc_p0_5(tc, ts, tw, i, m + 1) > eps)){
			m++;	
			p0 = calc_p0_5(tc, ts, tw, i, m);
		}
		printf("%d %d\n", i, m);
		pnm = pow(a, i + m) / (pow(i, m) * fact(i)) * p0;

		printf("pnm == %.10lf\n", pnm);

		K = a * (1.0 - pnm);
		L = ((pow(a, i + 1)*p0) / (i*fact(i))*((1.0 - pow(a / i, m)*(m + 1.0 - m*a / i)) / pow(1.0 - a / i, 2)));
		T = (L) * tc; 

		fprintf(k, "%d %lf\n",i, K);
		fprintf(l, "%d %lf\n",i, L);
		fprintf(t, "%d %lf\n",i, T);

		getchar();
		m = 1;
	}
}


void task_2(){
	double tc = 372.0;
	double ts = 12.0;
	double a = ts / tc;
	int N = 41;
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
		printf("task 2 : m == %d  p0 == %.10lf\n", m, p0);


		for(int k = 0; k <= N; k++){
			N_sr += (N - k) * PI[m][k] * p0;
		}
		printf("N_sr == %.20lf\n", N - N_sr);

		printf("pn == %.20lf\n", PI[m][N]);

		for(int i = 0; i <= m; i++){
			K_load += i * p0 * PI[m][i];
		}

		for(int j = m + 1; j <= N; j++){
			K_load += m * p0 * PI[m][j];
		}

		printf("K_load == %.20lf\n", K_load);
		K_load = 0.0;
		N_sr = 0.0;
		getchar();
	}
}

int main(int argc, char* argv[]){
	
	double tc, ts, tw; 
	int num_oper;
	FILE* file_oper_refuse;
	FILE* file_load_koef;
	FILE* file_reverse;
	FILE* file_load_koef2;
	FILE* file_len_queue;
	FILE* file_time_queue;
	FILE* file_load_koef5;
	FILE* file_len_queue5;
	FILE* file_time_queue5;


	if(argc < 4){
		fprintf(stderr, "Invalid arg num. Usage: tc, ts, tw\n");
		exit(INVALID_ARG);
	}

	tc = atoi(argv[1]);
	ts = atoi(argv[2]);
	tw = atoi(argv[3]);

	if((file_oper_refuse = fopen("file_oper_refuse.txt", "w")) == NULL){
		fprintf(stderr, "can't open file\n");
		exit(ERROR_OPEN_FILE);
	}

	if((file_load_koef = fopen("file_load_koef.txt", "w")) == NULL){
		fprintf(stderr, "can't open file\n");
		exit(ERROR_OPEN_FILE);
	}

	if((file_reverse = fopen("file_reverse.txt", "w")) == NULL){
		fprintf(stderr, "can't open file\n");
		exit(ERROR_OPEN_FILE);
	}

	if((file_load_koef2 = fopen("file_load_koef2.txt", "w")) == NULL){
		fprintf(stderr, "can't open file\n");
		exit(ERROR_OPEN_FILE);
	}

	if((file_len_queue = fopen("file_len_queue.txt", "w")) == NULL){
		fprintf(stderr, "can't open file\n");
		exit(ERROR_OPEN_FILE);
	}

	if((file_time_queue = fopen("file_time_queue.txt", "w")) == NULL){
		fprintf(stderr, "can't open file\n");
		exit(ERROR_OPEN_FILE);
	}

	if((file_load_koef5 = fopen("file_load_koef5.txt", "w")) == NULL){
		fprintf(stderr, "can't open file\n");
		exit(ERROR_OPEN_FILE);
	}

	if((file_len_queue5 = fopen("file_len_queue5.txt", "w")) == NULL){
		fprintf(stderr, "can't open file\n");
		exit(ERROR_OPEN_FILE);
	}

	if((file_time_queue5 = fopen("file_time_queue5.txt", "w")) == NULL){
		fprintf(stderr, "can't open file\n");
		exit(ERROR_OPEN_FILE);
	}


	assert(tc > 0 && ts > 0 && tw > 0 && "args must be positive\n");

	//num_oper = operators_num(tc, ts, 0.01, file_oper_refuse);

	//printf("%d\n", num_oper);

	//plot_load(tc, ts, num_oper, file_load_koef);
	//plot_reverse(tc, ts, num_oper, file_reverse, file_load_koef2);

	//find_serv_queue(tc, ts);

	//plot_len_time(tc, ts, num_oper, file_len_queue, file_time_queue);

	//printf("%lf\n", fact_rev(12, 12));

	//task_exit(tc, ts, tw, num_oper, file_load_koef5, file_len_queue5, file_time_queue5);

	task_2();

	fclose(file_oper_refuse);
	fclose(file_load_koef);
	return 0;
}


/*

Задача №1. Проектирование Call-центра.

Известно, что среднее время между звонками клиентов составляет Tc секунд, а среднее время обслуживания Ts секунд. Все потоки случайных событий считать пуассоновскими. Если все операторы заняты, звонок теряется.

1. Рассчитать минимальное число операторов, при котором доля отказов не превышает 25%, 10%, 5%, 1%. Построить график вероятности отказа от числа операторов. Построить график, иллюстрирующий коэффициент загрузки операторов в зависимости от их числа.

2. Построить график изменения доли отказов при замене каналов обслуживания местами в очереди (начиная с числа операторов, соответствующего 1% отказов в системе без очереди). Построить график, иллюстрирующий коэффициент загрузки операторов. Построить график математического ожидания длины очереди. Построить график, иллюстрирующий коэффициент занятости мест в очереди. Построить график математического ожидания времени пребывания клиентов в очереди.

3. Рассмотреть систему с ограниченной очередью. Подобрать такие параметры системы (число каналов обслуживания и число мест в очереди) чтобы коэффициенты загрузки операторов и занятости мест в очереди были максимально близки друг к другу.

4. Рассмотреть систему без ограничений на длину очереди. Построить график математического ожидания длины очереди в зависимости от числа операторов (вплоть до числа каналов, соответствующего 1% отказов в системе без очереди). Построить график, иллюстрирующий коэффициент загрузки операторов. Построить график математического ожидания времени пребывания клиентов в очереди.

5. Рассмотреть систему без ограничений на длину очереди, учитывающей фактор ухода клиентов из очереди (среднее приемлемое время ожидания – Tw секунд). Построить график средней длины очереди в зависимости от числа операторов (вплоть до числа каналов, соответствующего 1% отказов в системе без очереди). Построить график, иллюстрирующий коэффициент загрузки операторов. Построить график математического ожидания времени пребывания клиентов в очереди.

Задача №2. Проектирование производственного участка.

Имеется участок с N станками. Среднее время между наладками составляет Tc минут, среднее время наладки – Ts минут. Все потоки случайных событий считать пуассоновскими. Построить график зависимости числа простаивающих станков от числа наладчиков. Построить график зависимости коэффициента занятости наладчиков от их числа.

*/
