#ifndef __FACTORIAL__
#define __FACTORIAL__
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
#endif
