#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

typedef struct {
	double posX,posY,mass,velX,velY;
} body_t;

typedef struct BH {
	double corners[4]; // lower left corner, upper right corner
	struct BH **quadrants;
  body_t *body;

} BH_t;

void print_corners(BH_t *domain);

BH_t *create_domain(double llcx,
										double llcy,
										double urcx,
										double urcy)
{
	printf("create_domain params: %lf, %lf, %lf, %lf\n", llcx,llcy,urcx,urcy);
	BH_t *domain = malloc(sizeof(BH_t));
	domain->corners[0] = llcx;
	domain->corners[1] = llcy;
	domain->corners[2] = urcx;
	domain->corners[3] = urcy;
	domain->quadrants = malloc(4*sizeof(BH_t*));
	return domain;
}

void divide_domain(BH_t *domain) {
	domain->quadrants[0] = create_domain((domain->corners[0]+domain->corners[2])/2,
																			 (domain->corners[1]+domain->corners[3])/2,
																			 domain->corners[2],
																			 domain->corners[3]);

	domain->quadrants[1] = create_domain((domain->corners[0]+domain->corners[2])/2,
																			 domain->corners[1],
																			 domain->corners[2],
																			 (domain->corners[1]+domain->corners[3])/2);

	domain->quadrants[2] = create_domain(domain->corners[0],
																			 domain->corners[1],
																			 (domain->corners[0]+domain->corners[2])/2,
																			 (domain->corners[1]+domain->corners[3])/2);

	domain->quadrants[3] = create_domain(domain->corners[0],
																			 (domain->corners[1]+domain->corners[3])/2,
																			 (domain->corners[0]+domain->corners[2])/2,
																			 domain->corners[3]);
}

void print_corners(BH_t *domain) {
	printf("(%lf, %lf), (%lf, %lf)\n",
		domain->corners[0],
		domain->corners[1],
		domain->corners[2],
		domain->corners[3]);
}

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
	else
	{
		printf("which_quadrant: Something is very wrong.\n");
	}
	return quadrant;
}

void insert(body_t *body, BH_t *domain) {
	sleep(1);

	if (!domain->body) {
		domain->body = body;
	}
	else {

		divide_domain(domain);
		int q = which_quadrant(body,domain);
		printf("inserting new into Q%d:\n", q);
		insert(body, domain->quadrants[q]);
		
		q = which_quadrant(domain->body,domain);
		printf("inserting occupant into Q%d:\n", q);
		insert(domain->body, domain->quadrants[q]);
	}
}
/*
void main() {

	BH_t *nightsky = create_domain(0,0,1,1);
	body_t b1 = {0.375788, 0.506437,1,0,0};
	body_t b2 = {0.452506, 0.501224,1,0,0};
	insert(&b1,nightsky);
	insert(&b2,nightsky);
}
*/