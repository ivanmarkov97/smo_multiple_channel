#include "SMO.h"

#ifndef __SMO_LIMITED__
#define __SMO_LIMITED__

class SMO_LIMITED : public SMO {
public:
	SMO_LIMITED(double _tc, double _ts, double _tw): SMO(_tc, _ts, _tw) {};
	~SMO_LIMITED();

	void plot_len_time_work(int a, FILE* f0, FILE* f1, FILE* f2, FILE* f3, FILE* f4);
	void find_serv_queue();
};

#endif
