#include "SMO.h"

#ifndef __SMO_UNLIMITED__
#define __SMO_UNLIMITED__

class SMO_UNLIMITED : public SMO {
public:
	SMO_UNLIMITED(double _tc, double _ts, double _tw): SMO(_tc, _ts, _tw) {};
	~SMO_UNLIMITED();

	void plot_len_time(int, FILE*, FILE*, FILE*);
	void task_exit(int, FILE*, FILE*, FILE*);
	double calc_p0_5(int n, int m);
};

#endif
