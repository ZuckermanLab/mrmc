#ifndef SEDDD_H_INCLUDED
#define SEDDD_H_INCLUDED
//#define SEDDD
#include "ffield.h"
#ifdef SEDDD
//seddd_params has been moved to ffield.h.

void read_solvation_params(char * line, seddd_params * params);
#endif
#endif // SEDDD_H_INCLUDED
