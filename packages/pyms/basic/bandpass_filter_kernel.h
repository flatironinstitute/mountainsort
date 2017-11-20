/*
 * Copyright 2016-2017 Flatiron Institute, Simons Foundation
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef bandpass_filter_kernel_h
#define bandpass_filter_kernel_h

#include <stdint.h>

typedef int64_t bigint;

void bandpass_filter_kernel(bigint M,bigint N,float *X, double samplerate, double freq_min, double freq_max, double freq_wid, bigint start_write, bigint end_write);
void bandpass_filter_kernel_multithread(bigint M,bigint N,float *X, double samplerate, double freq_min, double freq_max, double freq_wid);

#endif