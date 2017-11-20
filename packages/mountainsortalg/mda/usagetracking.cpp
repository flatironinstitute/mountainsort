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
#include "usagetracking.h"
#include <QDebug>

static bigint num_files_open = 0;
static int64_t num_bytes_allocated = 0;
static bigint malloc_count = 0;
static bigint num_bytes_read = 0;
static bigint num_bytes_written = 0;

FILE* jfopen(const char* path, const char* mode)
{
    //printf("jfopen.\n");
    FILE* F = fopen(path, mode);
    if (!F)
        return 0;
    num_files_open++;
    return F;
}

void jfclose(FILE* F)
{
    //printf("jfclose.\n");
    if (!F)
        return;
    fclose(F);
    num_files_open--;
}

bigint jfread(void* data, size_t sz, bigint num, FILE* F)
{
    bigint ret = fread(data, sz, num, F);
    num_bytes_read += ret;
    return ret;
}

bigint jfwrite(void* data, size_t sz, bigint num, FILE* F)
{
    bigint ret = fwrite(data, sz, num, F);
    num_bytes_written += ret;
    return ret;
}

bigint jnumfilesopen()
{
    return num_files_open;
}

void* jmalloc(size_t num_bytes)
{
    //printf("jmalloc %d.\n",(bigint)num_bytes);
    if (!num_bytes)
        return 0;
    void* ret = malloc(num_bytes + 8); //add some space to track the number of bytes
    int32_t* tmp = (int32_t*)ret;
    tmp[0] = num_bytes;
    num_bytes_allocated += num_bytes;
    malloc_count++;
    return (void*)(((unsigned char*)ret) + 8);
}

void jfree(void* ptr)
{
    //printf("jfree.\n");
    if (!ptr)
        return;
    int64_t* tmp = (int64_t*)((unsigned char*)ptr - 8);
    int64_t num_bytes = tmp[0];
    free((void*)tmp);
    num_bytes_allocated -= num_bytes;
    malloc_count--;
}
bigint jmalloccount()
{
    return malloc_count;
}

int64_t jbytesallocated()
{
    return num_bytes_allocated;
}

bigint jnumbytesread()
{
    return num_bytes_read;
}

bigint jnumbyteswritten()
{
    return num_bytes_written;
}
