#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef __SMO__
#define __SMO__

class SMO {

protected:
	double tc;
	double ts;
	double tw;
	double a;
	double b;

public:
	SMO(double _tc = 0.0,
		double _ts = 0.0, 
		double _tw = 0.0): tc(_tc), ts(_ts), tw(_tw), a(_ts / _tc), b(_ts / _tw) {};
	~SMO();
	
	int num_oper(double, FILE*);
	void plot_load(int, FILE*);
	//void plot_len_time_work(int a, FILE* f0, FILE* f1, FILE* f2, FILE* f3, FILE* f4);
	void task_2(FILE*, FILE*);
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
	//void task_2();
};

#endif
