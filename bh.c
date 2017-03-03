#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
//#include "bh.h"

typedef struct {
  double posX,posY,mass,velX,velY;
} body_t;

typedef struct BH {
  double corners[4];
  struct BH **quadrants;
  body_t* body;
  double mass;
  double* com;
} BH_t;

void print_body(body_t* body);

void print_corners(BH_t* domain);

void print_domain_info(BH_t* domain);

void center_of_mass(double *b, double *c);

BH_t* create_domain(double llcx, double llcy, double urcx, double urcy) {

  BH_t *domain = malloc(sizeof(BH_t));
  domain->corners[0] = llcx;
  domain->corners[1] = llcy;
  domain->corners[2] = urcx;
  domain->corners[3] = urcy;
  domain->quadrants = malloc(4*sizeof(BH_t*));
  domain->body = NULL;
  domain->mass = 0;
  domain->com = malloc(2*sizeof(double));
  return domain;
}

void divide_domain(BH_t *domain) {
  double x[3],y[3];
  x[0] = domain->corners[0];
  x[2] = domain->corners[2];
  x[1] = (x[0]+x[2])/2;

  y[0] = domain->corners[1];
  y[2] = domain->corners[3];
  y[1] = (y[0]+y[2])/2;

  domain->quadrants[0] = create_domain(x[1],y[1],x[2],y[2]);
  domain->quadrants[1] = create_domain(x[1],y[0],x[2],y[1]);
  domain->quadrants[2] = create_domain(x[0],y[0],x[1],y[1]);
  domain->quadrants[3] = create_domain(x[0],y[1],x[1],y[2]);
}

/* Determine which quadrant of domain that body belongs to */
int which_quadrant(body_t *body, BH_t *domain) {

  int quadrant;
  if (body->posX <= domain->corners[2]
&& body->posY <= domain->corners[3]
&& body->posX > (domain->corners[0]+domain->corners[2])/2
&& body->posY > (domain->corners[1]+domain->corners[3])/2)
  {
    quadrant = 0;
  }
  else if (body->posX <= domain->corners[2]
&& body->posY <= (domain->corners[1]+domain->corners[3])/2
&& body->posX > (domain->corners[0]+domain->corners[2])/2
&& body->posY > domain->corners[1])
  {
    quadrant = 1;
  }
  else if (body->posX <= (domain->corners[0]+domain->corners[2])/2
&& body->posY <= (domain->corners[1]+domain->corners[3])/2
&& body->posX > domain->corners[0]
&& body->posY > domain->corners[1])
  {
    quadrant = 2;
  }
  else if (body->posX <= (domain->corners[0]+domain->corners[2])/2
&& body->posY <= domain->corners[3]
&& body->posX > domain->corners[0]
&& body->posY > (domain->corners[1]+domain->corners[3])/2)
  {
    quadrant = 3;
  }
  else {
    printf("Body does not belong to domain.\n");
  }

  return quadrant;
}

void insert(body_t *body, BH_t *domain) {
  sleep(1);

  if (!domain->body) {
    domain->body = body;
    domain->mass += body->mass;
    printf("Body inserted:\t"); print_body(body);
    printf("In domain:\t"); print_corners(domain);
    printf("This domain now has mass %lf\n", domain->mass);
    // update mass and com?
  }
  else {
    /* Divide the domain, put body in its appropriate quadrant. */
    divide_domain(domain);
    // update mass and com?
    int q = which_quadrant(body, domain);
    insert(body, domain->quadrants[q]); // crash happens on this line
    domain->mass += body->mass;
    printf("Super-domain now has mass: %lf\n", domain->mass);
    printf("--------------------------------------\n");

    /* Move existing body to its appropriate quadrant,
       and set domain's body to null. */
    q = which_quadrant(domain->body, domain);
    insert(domain->body, domain->quadrants[q]);
//    domain->body = NULL;
  }
}

double distance(body_t *b, body_t *c) {

  double dx = b->posX - c->posX;
  double dy = b->posY - c->posY;
  return sqrt(dx*dx+dy*dy);
}

void print_body(body_t* body) {
  printf("(%lf,%lf)\n", body->posX, body->posY);
}

void print_corners(BH_t* domain) {
  printf("%lf,%lf,%lf,%lf\n", domain->corners[0],domain->corners[1],
domain->corners[2],domain->corners[3]);
}

void print_domain_info(BH_t* domain) {
  printf("Corners: "); print_corners(domain);

  printf("Has quadrants? ");
  if (!domain->quadrants[0])
    printf("No\n");
  else
    printf("Yes\n");

  printf("Has body? ");
  if (!domain->body)
    printf("No\n");
  else
    printf("Yes\n");

  printf("Mass: %lf\n", domain->mass);
}

int main(int argc, char** argv) {

  // Create root domain D:
  BH_t* domain = create_domain(0,0,1,1);

  // Root domain now has two null pointers:
//  printf("Addr to quadrants array: %p\n", domain->quadrants);
//  printf("Addr to body: %p\n", domain->body);

  // Create 4 bodies, one for each quadrant:
  body_t b1 = {0.87,0.87,1,0,0}; // D->Q1
  body_t b2 = {0.87,0.37,1,0,0}; // D->Q2
  body_t b3 = {0.37,0.37,1,0,0}; // D->Q3
  body_t b4 = {0.37,0.87,1,0,0}; // D->Q4

  // Insert the first body:
  printf("Inserting b1:\n");
  insert(&b1, domain);
  //print_domain_info(domain);

  // At this point, domain has a body, but still no quadrants:
//  printf("Addr to quadrants array: %p\n", domain->quadrants);

  // Insert another body:
  //insert(&b2, domain);

  /* This caused a segfault. Thus, we need to instantiate
     the quadrants array from the start. malloc in create_domain.*/

  // Now run the code again. We now have a segfault in insert.

  // Comment out the second insert above.

  // Divide the domain manually:
  //divide_domain(domain);

  // Attempt access to domain->body:
//  if (!domain->body)
//    printf("success\n");

  // So far so good. Comment out the manual divide.

  // Try the second insertion again:
//  insert(&b2, domain);

  // Still segfault. Comment out 2nd insert.

  // Manually divide the domain, and then attempt to access 2nd quadrant:
//  divide_domain(domain);
//  if (domain->quadrants[1])
//    printf("success\n");

  // This worked. Comment it out.

  // Problem was in which_quadrant: it didn't return anything. Now fixed.

  // Now retry adding the 2nd body:
//  insert(&b2, domain);

  // Still segfault.

  // Added some prints and sleep. Retry:
//  insert(&b2, domain);

  // There was an error in insert, now fixed.

  // Try inserting the remaining bodies:
  printf("Inserting b2:\n");
  insert(&b2, domain);
//  print_domain_info(domain);
//  print_domain_info(domain->quadrants[1]);
  printf("Inserting b3:\n");
  insert(&b3, domain);
//  print_domain_info(domain);
//  print_domain_info(domain->quadrants[2]);
  printf("Inserting b4:\n");
  insert(&b4, domain);
//  print_domain_info(domain);
//  print_domain_info(domain->quadrants[3]);

  // So far so good. Added mass accumulator.

  // Now define a point that will go in root domain's 1st quadrant:
  body_t b11 = {0.87,0.63,1,0,0}; // D->Q1->Q2

  // Insert it:
  insert(&b11, domain);
  print_domain_info(domain);
  print_domain_info(domain->quadrants[0]);

  // That went well. Now print info about all domains:
/*  printf("Domain:\n"); print_domain_info(domain);
  printf("\n");
  printf("D->Q1\n"); print_domain_info(domain->quadrants[0]);
  printf("\n");
  printf("D->Q2\n"); print_domain_info(domain->quadrants[1]);
  printf("\n");
  printf("D->Q3\n"); print_domain_info(domain->quadrants[2]);
  printf("\n");
  printf("D->Q4\n"); print_domain_info(domain->quadrants[3]);
*/
  /* There is at least one error: reported mass.
     All of them should have mass at least 1. */

  // Print info after insertions.

}












