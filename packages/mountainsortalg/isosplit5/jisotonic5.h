/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
*******************************************************/

#ifndef jisotonic5_h
#define jisotonic5_h

#include <stdlib.h>

typedef int64_t bigint;

void jisotonic5(bigint N, float* BB, float* MSE, float* AA, float* WW);
void jisotonic5_updown(bigint N, float* out, float* in, float* weights);
void jisotonic5_downup(bigint N, float* out, float* in, float* weights);
void jisotonic5_sort(bigint N, float* out, const float* in);

#endif
