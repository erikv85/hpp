#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPS 1e-03
#define DATA_LEN

void test() {

	printf("include-test: %lf\n", pow(2,3));
}

void calculate_force(double *data, int nparticles, int body, double *force) {

	int k;
	double rx,ry,r,denom,factor;
	for (k=0; k<nparticles*DATA_LEN; k+=DATA_LEN) {

		rx = data[body]-data[k];
		ry = data[body+1]-data[k+1];
		r = sqrt( rx*rx + ry*ry );
		denom = cube(r,EPS);
		factor = data[k+2]/denom;
		force[0] += rx*factor;
		force[1] += ry*factor;
	}
}

void move(double *data, double *buffer, double *force, int body, double dt) {
	buffer[body+3] = data[body+3]+dt*force[0];
	buffer[body+1+3] = data[body+1+3]+dt*force[1];
	buffer[body+2] = data[body+2];
	buffer[body] = data[body]+dt*buffer[body+3];
	buffer[body+1] = data[body+1]+dt*buffer[body+1+3];
}