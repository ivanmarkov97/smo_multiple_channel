#include "SMO.h"

#ifndef __SMO_LIMITED__
#define __SMO_LIMITED__

class SMO_LIMITED : public SMO {
public:
	SMO_LIMITED(double _tc, double _ts, double _tw): SMO(_tc, _ts, _tw) {};
	~SMO_LIMITED();

	void find_serv_queue();
};

#endif
