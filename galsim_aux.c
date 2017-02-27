#include <time.h>

struct timespec {
	time_t tv_sec;
	long tv_nsec;
};

int end_of_first_half(int n) {

  if (2*(n/2)!=n)
    return -1;
  else
    return n/2;
}

timespec diff(timespec start, timespec end) {

	timespec temp;
	if ((end.tv_nsec-start.tv_nsec)<0) {
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	}
	else {
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp;	
}