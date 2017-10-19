#ifndef bandpass_filter_kernel_h
#define bandpass_filter_kernel_h

#include <stdint.h>

typedef int64_t bigint;

void bandpass_filter_kernel(bigint M,bigint N,float *X, double samplerate, double freq_min, double freq_max, double freq_wid, bigint start_write, bigint end_write);
void bandpass_filter_kernel_multithread(bigint M,bigint N,float *X, double samplerate, double freq_min, double freq_max, double freq_wid);

#endif