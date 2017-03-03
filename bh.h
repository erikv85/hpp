// STRUCTS
typedef struct {
  double posX,posY,mass,velX,velY;
} body_t;

typedef struct BH {
  double corners[4];
  struct BH** quadrants;
  body_t* body;
  double mass;
  double* com;
} BH_t;

// FUNCTION HEADERS
void print_body(body_t* body);

void print_corners(BH_t* domain);

void print_domain_info(BH_t* domain);

void center_of_mass(double* b, double* c);

BH_t* create_domain(double llcx, double llcy, double urcx, double urcy);

void divide_domain(BH_t* domain);

int which_quadrant(body_t* body, BH_t* domain);

void insert(body_t* body, BH_t* domain);

double distance(body_t* b, body_t* c);
