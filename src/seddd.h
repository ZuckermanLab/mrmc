#ifndef SEDDD_H_INCLUDED
#define SEDDD_H_INCLUDED

#include "ffield.h"
#define SEDDD

struct seddd_params {
    double c, eps0, eps1;
    double hydration_volume[MAX_NUM_OF_ATOM_CLASSES];
    double hydration_shell_thickness[MAX_NUM_OF_ATOM_CLASSES];
};

void read_solvation_params(char * line, seddd_params * params);
#endif // SEDDD_H_INCLUDED
