/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
*******************************************************/

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
