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
	int num_oper;
	FILE* file_oper_refuse;
	FILE* file_load_koef;
	FILE* file_refuse_reverse;
	FILE* file_load_koef2;
	FILE* file_time_queue2;
	FILE* file_len_queue2;
	FILE* file_work_queue2;
	FILE* file_load_koef4;
	FILE* file_len_queue;
	FILE* file_time_queue;
	FILE* file_load_koef5;
	FILE* file_len_queue5;
	FILE* file_time_queue5;
	FILE* file_free_works;
	FILE* file_work_koef;

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

	if((file_refuse_reverse = fopen("file_refuse_reverse.txt", "w")) == NULL){
		fprintf(stderr, "can't open file\n");
		exit(ERROR_OPEN_FILE);
	}

	if((file_load_koef2 = fopen("file_load_koef2.txt", "w")) == NULL){
		fprintf(stderr, "can't open file\n");
		exit(ERROR_OPEN_FILE);
	}

	if((file_len_queue2 = fopen("file_len_queue2.txt", "w")) == NULL){
		fprintf(stderr, "can't open file\n");
		exit(ERROR_OPEN_FILE);
	}

	if((file_time_queue2 = fopen("file_time_queue2.txt", "w")) == NULL){
		fprintf(stderr, "can't open file\n");
		exit(ERROR_OPEN_FILE);
	}

	if((file_work_queue2 = fopen("file_work_queue2.txt", "w")) == NULL){
		fprintf(stderr, "can't open file\n");
		exit(ERROR_OPEN_FILE);
	}

	if((file_load_koef4 = fopen("file_load_koef4.txt", "w")) == NULL){
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

	if((file_free_works = fopen("file_free_works.txt", "w")) == NULL){
		fprintf(stderr, "can't open file\n");
		exit(ERROR_OPEN_FILE);
	}

	if((file_work_koef = fopen("file_work_koef.txt", "w")) == NULL){
		fprintf(stderr, "can't open file\n");
		exit(ERROR_OPEN_FILE);
	}

	assert(tc > 0 && ts > 0 && tw > 0 && "args must be positive\n");

	SMO *smo = new SMO(tc, ts, tw);
	SMO_LIMITED *smo_lim = new SMO_LIMITED(tc, ts, tw);
	SMO_UNLIMITED *smo_unlim = new SMO_UNLIMITED(tc, ts, tw);

	int num = smo->num_oper(0.01, file_oper_refuse);
	smo->plot_load(num, file_load_koef);
	smo->plot_len_time_work(num, file_refuse_reverse, file_load_koef2, file_len_queue2, file_work_queue2, file_time_queue2);
	smo_lim->find_serv_queue();
	smo_unlim->plot_len_time(num, file_load_koef4, file_len_queue, file_time_queue);
	smo_unlim->task_exit(num, file_load_koef5, file_len_queue5, file_time_queue5);
	smo->task_2(file_free_works, file_work_koef);

	return 0;
}