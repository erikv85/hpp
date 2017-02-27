#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
//#include "graphics/graphics.c"
#include "file_operations/file_operations.c"
//#include "galsim_aux.c"
//#include "physics.c"
//#include "bh.c"

#define EPS 1e-03
#define DATA_LEN 5

double cube(double x, double y);

void calculate_force(double *data, int nparticles, int body, double *force);

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
  double *force = malloc(2*sizeof(double));
  double *buffer = malloc(N*DATA_LEN*sizeof(double));
  double *tmp;

  // Load data:
  double *data = malloc(N*DATA_LEN*sizeof(double));
  if (read_doubles_from_file(N*DATA_LEN, data, filename)<0) {
    perror(NULL);
    return -1;
  }

/*
  // Format data and create Barnes-Hut domain:
  body_t **fmtd_data = malloc(N*sizeof(body_t*));
  BH_t *domain = create_domain(0,0,1,1);
  int i;
  for (i=0; i<N; i++) {
    fmtd_data[i] = format_data(&data[i*DATA_LEN]);
    printf("inserting data i: %d\n", i);
    printf("body coordinates: %0.20lf, %0.20lf\n", fmtd_data[i]->posX, fmtd_data[i]->posY);
    insert(fmtd_data[i],domain);
  }
*/

  clock_t outer_start,outer_end;
  int i,j,k;
  /* WORK */
  outer_start = clock();
  for (i=0; i<nsteps; i++) {

    for (j=0; j<N*DATA_LEN; j+=DATA_LEN) {

      force[0]=0;
      force[1]=0;
      
      calculate_force(data,N,j,force);
      //printf("j, forceX: %d, %lf\n", j,force[j]);
      force[0] = -G*force[0];
      force[1] = -G*force[1];
      //printf("forceX: %lf\n", force[0]);

      move(data,buffer,force,j,dt);
      printf("j, data[0]: %d, %lf\n", j, data[0]);
      
    }
    //printf("data[0]: %lf\n", data[0]);
    //sleep(1);

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
  free(force);
  free(buffer);

  return 0;
}

double cube(double x, double y) {

  return x*x*x+3*x*x*y+3*x*y*y+y*y*y;
}

void calculate_force(double *data, int nparticles, int body, double *force) {

  int k;
  double rx,ry,r,denom,factor;
  for (k=0; k<nparticles*DATA_LEN; k+=DATA_LEN) {
    //printf("k, data[k]: %d, %lf\n", k,data[k]);
    rx = data[body]-data[k];  //printf("rx: %lf\n", rx);
    ry = data[body+1]-data[k+1];
    r = sqrt( rx*rx + ry*ry );
    denom = cube(r,EPS);
    factor = data[k+2]/denom;
    force[0] += rx*factor;
    force[1] += ry*factor;
    //printf("k, forceX: %d, %lf\n", k,force[0]);
    //sleep(0.5);
  }
}

void move(double *data, double *buffer, double *force, int body, double dt) {
  buffer[body+3] = data[body+3]+dt*force[0];
  buffer[body+1+3] = data[body+1+3]+dt*force[1];
  buffer[body+2] = data[body+2];
  buffer[body] = data[body]+dt*buffer[body+3];
  buffer[body+1] = data[body+1]+dt*buffer[body+1+3];
}