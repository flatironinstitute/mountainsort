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

#ifndef USAGETRACKING
#define USAGETRACKING

#include <stdlib.h>
#include <stdio.h>
typedef int64_t bigint;

FILE* jfopen(const char* path, const char* mode);
void jfclose(FILE* F);
bigint jfread(void* data, size_t sz, bigint num, FILE* F);
bigint jfwrite(void* data, size_t sz, bigint num, FILE* F);
bigint jnumfilesopen();

void* jmalloc(size_t num_bytes);
void jfree(void* ptr);
bigint jmalloccount();
int64_t jbytesallocated();
bigint jnumbytesread();
bigint jnumbyteswritten();

#endif // USAGETRACKING
