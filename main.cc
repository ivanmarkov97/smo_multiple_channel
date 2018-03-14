#include "SMO.h"
#include "SMO_LIMITED.h"
#include "SMO_UNLIMITED.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define INVALID_ARG -1.0
#define ERROR_OPEN_FILE -2.0

int main(int argc, char* argv[]){
	
	double tc, ts, tw; 

	FILE* files[15];
	const char* file_names[15] = {
		"file_oper_refuse.txt",
        "file_load_koef.txt",
        "file_refuse_reverse.txt",
        "file_load_koef2.txt",
        "file_len_queue2.txt",
        "file_work_queue2.txt",
        "file_time_queue2.txt",
        "file_load_koef4.txt",
        "file_len_queue.txt",
        "file_time_queue.txt",
        "file_load_koef5.txt",
        "file_len_queue5.txt",
        "file_time_queue5.txt",
        "file_free_works.txt",
        "file_work_koef.txt"   
	};

	if(argc < 4){
		fprintf(stderr, "Invalid arg num. Usage: tc, ts, tw\n");
		exit(INVALID_ARG);
	}

	tc = atoi(argv[1]);
	ts = atoi(argv[2]);
	tw = atoi(argv[3]);

	for(int i = 0; i < 15; i++){
		if((files[i] = fopen(file_names[i], "w")) == NULL){
			fprintf(stderr, "can't open file\n");
			exit(ERROR_OPEN_FILE);
		}
	}

	assert(tc > 0 && ts > 0 && tw > 0 && "args must be positive\n");

	SMO *smo = new SMO(tc, ts, tw);
	SMO_LIMITED *smo_lim = new SMO_LIMITED(tc, ts, tw);
	SMO_UNLIMITED *smo_unlim = new SMO_UNLIMITED(tc, ts, tw);

	int num = smo->num_oper(0.01, files[0]);
	smo->plot_load(num, files[1]);
	smo_lim->plot_len_time_work(num, files[2], files[3], files[4], files[5], files[6]);
	smo_lim->find_serv_queue();
	smo_unlim->plot_len_time(num, files[7], files[8], files[9]);
	smo_unlim->task_exit(num, files[10], files[11], files[12]);
	smo->task_2(files[13], files[14]);

	return 0;
}