#include <math.h>
#include <stdio.h>
/*
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.keio.ac.jp/matumoto/emt.html
   email: matumoto@math.keio.ac.jp
*/

#include "mt.h"

// initializes mt[NN] with a seed //
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<NN; mti++) {
        mt[mti] =
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        // See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier.
        // In the previous versions, MSBs of the seed affect
        // only MSBs of the array mt[].
        mt[mti] &= 0xffffffffUL;
        // for >32 bit machines
    }
}

// generates a random number on [0,0xffffffff]-interval
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    // mag01[x] = x * MATRIX_A  for x=0,1

    if (mti >= NN) { // generate NN words at one time
        int kk;

        if (mti == NN+1)         // if init_genrand() has not been called,
	  init_genrand(5489UL); // a default initial seed is used

        for (kk=0;kk<NN-MM;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+MM] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<NN-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(MM-NN)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[NN-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[NN-1] = mt[MM-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    // Tempering
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
    //printf("Random number: %.8lx\n",y);
    return y;
}

// generates a random number on [0,1]-real-interval
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0);
    // divided by 2^32-1
}

// generates a random number on [0,1)-real-interval
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0);
    // divided by 2^32
}

// generates a random number on (0,1)-real-interval
double genrand_real3(void)
{
    double r;
    r=(((double)genrand_int32()) + 0.5)*(1.0/4294967296.0);
    return r;
    // divided by 2^32
}

// generates a random number on [0,1) with 53-bit resolution
double genrand_res53(void)
{
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6;
    return(a*67108864.0+b)*(1.0/9007199254740992.0);
}

// generates a gaussian random number using genrand_real3()
double gauss(double gw){
  int iset = 0;
  double gset;
  double fac,rsq,v1,v2;
  if(iset==0){
    do{
      v1 = 2.0*genrand_real3()-1.0; //from a uniform dist (-1,1)
      v2 = 2.0*genrand_real3()-1.0;
      rsq = v1*v1+v2*v2;
    }
    while(rsq>=1.0 || rsq==0.0);
    fac = sqrt(-2.0*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;
    return v2*fac*(gw);    // without () term width is -1 to 1
  }
  else{
      iset = 0;
      return gset*(gw);      // without () term width is -1 to 1
  }
}

//justin's code, these have to be inside mt.cpp in order to reference the mt[] array

void write_rng_state(FILE * output)
{
	int i;
	fprintf(output,"%d\n",mti);
	for (i=0; i<NN; i+=8) fprintf(output,"%.8lx %.8lx %.8lx %.8lx %.8lx %.8lx %.8lx %.8lx\n",mt[i],mt[i+1],mt[i+2],mt[i+3],mt[i+4],mt[i+5],mt[i+6],mt[i+7]);
}

void read_rng_state(FILE * input)
{
	int i;
	fscanf(input,"%d\n",&mti);
        for (i=0; i<NN; i+=8) fscanf(input,"%lx %lx %lx %lx %lx %lx %lx %lx\n",&mt[i],&mt[i+1],&mt[i+2],&mt[i+3],&mt[i+4],&mt[i+5],&mt[i+6],&mt[i+7]);
}	
