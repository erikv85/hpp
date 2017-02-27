#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../file_operations/file_operations.c" // changed

static void update_maxdiff(double dx, double dy, double* maxabsdiff) {
  double absdiff = sqrt(dx*dx+dy*dy);
  if(absdiff > *maxabsdiff)
    *maxabsdiff = absdiff;
}

/* The idea with the check_that_numbers_seem_OK() function is to check
   that there are no strange numbers like "nan" that may give problems
   when we try to compare the numbers later. */
int check_that_numbers_seem_OK(int n, double* buf) {
  const double minAllowedValue = -1e10;
  const double maxAllowedValue = 1e10;
  int OKcount = 0;
  int i;
  for(i = 0; i < n; i++) {
    double a = buf[i];
    if(a >= minAllowedValue && a <= maxAllowedValue)
      OKcount++;
  }
  if(OKcount == n)
    return 0;
  else
    return -1;
}

int main(int argc, const char* argv[]) {
  if(argc != 4) {
    printf("Give 3 input args: N gal1.gal gal2.gal\n");
    return -1;
  }
  int N = atoi(argv[1]);
  const char* fileName1 = argv[2];
  const char* fileName2 = argv[3];
  printf("N = %d\n", N);
  printf("fileName1 = '%s'\n", fileName1);
  printf("fileName2 = '%s'\n", fileName2);
  /* Read files. */
  double buf1[5*N];
  if(read_doubles_from_file(5*N, buf1, fileName1) != 0) {
    printf("Error reading file '%s'\n", fileName1);
    return -1;
  }
  if(check_that_numbers_seem_OK(5*N, buf1) != 0) {
    printf("Error: strange numbers found in file '%s'.\n", fileName1);
    return -1;
  }
  double buf2[5*N];
  if(read_doubles_from_file(5*N, buf2, fileName2) != 0) {
    printf("Error reading file '%s'\n", fileName2);
    return -1;
  }
  if(check_that_numbers_seem_OK(5*N, buf2) != 0) {
    printf("Error: strange numbers found in file '%s'.\n", fileName2);
    return -1;
  }
  /* Compare positions and velocities. */
  double pos_maxdiff = 0;
  double vel_maxdiff = 0;
  int i;
  for(i = 0; i < N; i++) {
    double pos_dx = buf1[i*5+0] - buf2[i*5+0];
    double pos_dy = buf1[i*5+1] - buf2[i*5+1];
    double vel_dx = buf1[i*5+3] - buf2[i*5+3];
    double vel_dy = buf1[i*5+4] - buf2[i*5+4];
    update_maxdiff(pos_dx, pos_dy, &pos_maxdiff);
    update_maxdiff(vel_dx, vel_dy, &vel_maxdiff);
  }
  printf("pos_maxdiff = %16.12f\n", pos_maxdiff);
  printf("vel_maxdiff = %16.12f\n", vel_maxdiff);
  return 0;
}
