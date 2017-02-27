#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "file_operations/file_operations.c"
#define EPS 1e-03
#define DATA_LEN 5

double cube(double x, double y);

void calculate_force(double *data,
                    int nparticles,
                    //int body,
                    double *force,
                    int indA,
                    int indB,
                    double G,
                    double *buffer,
                    double dt);

void move(double *data, double *buffer, double *force, int body, double dt);

int main(int argc, char *argv[]) {

  // Handle program arguments:
  int N, nsteps, graphics, nthreads;
  char *filename;
  double dt, theta_max;
  if (argc<7 || argc>8) {
    printf("Usage: N, filename, dt, theta_max,graphics\n");
    return -1;
  }
  N = atoi(argv[1]);
  filename = argv[2];
  nsteps = atoi(argv[3]);
  dt = atof(argv[4]);
  theta_max = atof(argv[5]);
  graphics = atoi(argv[6]);
  if (argc==8) {
    nthreads = atoi(argv[7]);
    printf("Barnes-Hut calculation with %d thread(s).\n", nthreads);
  }
  else {
    printf("Barnes-Hut calculation with 1 thread.\n");
  }

  // Physics parameters:
  double G = 100/((double)N);
  double *force = malloc(N*DATA_LEN*sizeof(double));
  double *buffer = malloc(N*DATA_LEN*sizeof(double));
  double *tmp;

  // Load data:
  double *data = malloc(N*DATA_LEN*sizeof(double));
  if (read_doubles_from_file(N*DATA_LEN, data, filename)<0) {
    perror(NULL);
    return -1;
  }

  nthreads = 2;
  int indA_1 = 0;
  int indB_1 = N/2;
  int indA_2 = N/2;
  int indB_2 = N;
  printf("Thread 1's indices: %d, %d\n", indA_1, indB_1);
  printf("Thread 2's indices: %d, %d\n", indA_2, indB_2);
  printf("\n\n");

  clock_t outer_start,outer_end;
  int i,j,k;
  /* WORK */
  outer_start = clock();
  for (i=0; i<nsteps; i++) { printf("----NSTEPS=%d-----\n", i);
/*
    for (j=indA_2; j<indB_2*DATA_LEN; j+=DATA_LEN) {

      force[j] = 0;
      force[j+1] = 0;
      double rx,ry,r,denom,factor;
      for (k=0; k<N*DATA_LEN; k+=DATA_LEN) {
        
        rx = data[j]-data[k];
        ry = data[j+1]-data[k+1];
        r = sqrt( rx*rx + ry*ry );
        denom = cube(r,EPS);
        factor = data[k+2]/denom;
        force[j] += rx*factor;
        force[j+1] += ry*factor;

      }
      force[j] = -G*force[j];
      force[j+1] = -G*force[j+1];

      move(data,buffer,force,j,dt);

    }
*/
    //calculate_force(data,N,force,0,N,G,buffer,dt);

    #pragma omp parallel num_threads(nthreads)
    {
      int id = omp_get_thread_num();
      if (id==0) {
		printf("thread %d will handle indices %d<=j<%d\n", omp_get_thread_num(), indA_1*DATA_LEN, indB_1*DATA_LEN);
        calculate_force(data,N,force,indA_1,indB_1,G,buffer,dt);
	  }
      else {
		//usleep(100);
		printf("thread %d will handle indices %d<=j<%d\n", omp_get_thread_num(), indA_2*DATA_LEN, indB_2*DATA_LEN);
        calculate_force(data,N,force,indA_2,indB_2,G,buffer,dt);
      }
    }

    tmp = data;
    data = buffer;
    buffer = tmp;
  }
  outer_end = clock();
  printf("Outer time (s):\t%lf\n", (outer_end-outer_start)/(double)CLOCKS_PER_SEC);
  
  // Save results to file:
  char saveto[] = "result.gal";
  if (write_doubles_to_file(N*DATA_LEN, data, saveto)<0) {
    perror(NULL);
    return -1;
  }

  // Clear heap:
  free(data);
  //free(force);
  //free(buffer);

  return 0;
}

double cube(double x, double y) {

  return x*x*x+3*x*x*y+3*x*y*y+y*y*y;
  printf("%lf\n", y);
}

void calculate_force(double *data,
                    int nparticles,
                    //int body,
                    double *force,
                    int indA,
                    int indB,
                    double G,
                    double *buffer,
                    double dt) {
 
  int j,k;
  double rx,ry,r,denom,factor;
  for (j=indA*DATA_LEN; j<indB*DATA_LEN; j+=DATA_LEN) {
	printf("thread %d got indices %d and %d\n", omp_get_thread_num(), indA, indB);
    force[j] = 0;
    force[j+1] = 0;
    for (k=0; k<nparticles*DATA_LEN; k+=DATA_LEN) {

      rx = data[j]-data[k];
      ry = data[j+1]-data[k+1];
      r = sqrt( rx*rx + ry*ry );
      denom = cube(r,EPS);
      factor = data[k+2]/denom;
      force[j] += rx*factor;
      force[j+1] += ry*factor;
    }
    force[j] = -G*force[j];
    force[j+1] = -G*force[j+1];
	printf("thread %d handling index %d (body %d), calling move:\n", omp_get_thread_num(), j, j/DATA_LEN);
    move(data,buffer,force,j,dt);
    //printf("j, data[0]: %d, %lf\n", j, data[0]);
  }
}

void move(double *data, double *buffer, double *force, int body, double dt) {
  buffer[body+3] = data[body+3]+dt*force[body];			// update velX
  buffer[body+1+3] = data[body+1+3]+dt*force[body+1];	// update velY
  buffer[body+2] = data[body+2];						// copy mass
  buffer[body] = data[body]+dt*buffer[body+3];			// update posX
  buffer[body+1] = data[body+1]+dt*buffer[body+1+3];	// update posY
  int x = omp_get_thread_num();
  printf("thread %d moved body %d to: %lf (posX)\n", omp_get_thread_num(), body/DATA_LEN, buffer[body]);
}
