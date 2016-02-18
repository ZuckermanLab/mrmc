//begin mersenne twister header

#ifndef __MT_H__
#define __MT_H__


//#include <stdio.h>

// Period parameters  
#define NN 624
#define MM 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[NN]; /* the array for the state vector  */
static int mti=NN+1; /* mti==NN+1 means mt[NN] is not initialized */

void init_genrand(unsigned long);
unsigned long genrand_int32(void);
double genrand_real1(void);
double genrand_real2(void);
double genrand_real3(void);
double genrand_res53(void);
//end mersenne twister header

//begin gauss header
//#include <math.h>
double gauss(double);
//end gauss header


//justin's code
void write_rng_state(FILE * output);
void read_rng_state(FILE * input);
#endif //__MT_H__
