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

#ifndef MDAIO_H
#define MDAIO_H

#include <stdlib.h>
typedef int64_t bigint;

//#include "mlcommon.h"

/*
 * For a discussion see: http://wisdmhub.org/dev2/apps/mdaformat/
 *
 * This is a C function.
 * To include within a C++ program use
 * #include "mdaio.h"
 * You will also need to link in mdaio.cpp
 *
 * Example usage: See the transpose_array function at the bottom of this file
 *	(and note that this would be more elegant with C++) ;)
 *
 * Another example usage: scda_ss/processing/filters/expfilter/expfilter.c
 *
 * The key trick here is the convenience of the functions
 * mda_read_float32, mda_read_byte, ..., and mda_write_float, mda_write_byte, ...
 * because they allows you to read/write of a chosen type without worrying about
 * the actual data storage type. Normally you would need to have if/else statements
 * to handle every possible type. Think about how horrible this would be for even
 * the simple transpose_array function below.
 *
 * Note that we don't handle the MDA_TYPE_COMPLEX case in this file
 *
 * */

#include <stdio.h>
#include <stdlib.h>

#define MDAIO_MAX_DIMS 50
#define MDAIO_MAX_SIZE 1e15
#define MDAIO_TYPE_COMPLEX -1
#define MDAIO_TYPE_BYTE -2
#define MDAIO_TYPE_FLOAT32 -3
#define MDAIO_TYPE_INT16 -4
#define MDAIO_TYPE_INT32 -5
#define MDAIO_TYPE_UINT16 -6
#define MDAIO_TYPE_FLOAT64 -7
#define MDAIO_TYPE_UINT32 -8
//FUTURE:
//#define MDAIO_TYPE_INT64 -9
//#define MDAIO_TYPE_UINT64 -10

#include <stdint.h>

struct MDAIO_HEADER {
    //note that all these data are *signed* 32-bit integers -- a historical choice
    int32_t data_type; //the data type in the file
    int32_t num_bytes_per_entry; //number of bytes per data entry (redundant)
    int32_t num_dims; //the number of dimensions to read/write
    uint64_t dims[MDAIO_MAX_DIMS]; //the size of the dimensions

    //the following make sense for reading but not writing -- purposes of info
    //note that the header size will be (3+num_dims)x4 or 3x4+num_dims*8
    bigint header_size; //the size of the header in bytes
};

//simply read, write or copy the mda header
bigint mda_read_header(struct MDAIO_HEADER* H, FILE* input_file);
bigint mda_write_header(struct MDAIO_HEADER* H, FILE* output_file);
void mda_copy_header(struct MDAIO_HEADER* Hdst, const struct MDAIO_HEADER* Hsrc);

//the following can be used no matter what the underlying data type is
bigint mda_read_byte(unsigned char* data, struct MDAIO_HEADER* H, bigint n, FILE* input_file);
bigint mda_read_float32(float* data, struct MDAIO_HEADER* H, bigint n, FILE* input_file);
bigint mda_read_int16(int16_t* data, struct MDAIO_HEADER* H, bigint n, FILE* input_file);
bigint mda_read_int32(int32_t* data, struct MDAIO_HEADER* H, bigint n, FILE* input_file);
bigint mda_read_uint16(uint16_t* data, struct MDAIO_HEADER* H, bigint n, FILE* input_file);
bigint mda_read_float64(double* data, struct MDAIO_HEADER* H, bigint n, FILE* input_file);
bigint mda_read_uint32(uint32_t* data, struct MDAIO_HEADER* H, bigint n, FILE* input_file);

//the following can be used no matter what the underlying data type is
bigint mda_write_byte(unsigned char* data, struct MDAIO_HEADER* H, bigint n, FILE* output_file);
bigint mda_write_float32(const float* data, struct MDAIO_HEADER* H, bigint n, FILE* output_file);
bigint mda_write_int16(int16_t* data, struct MDAIO_HEADER* H, bigint n, FILE* output_file);
bigint mda_write_int32(int32_t* data, struct MDAIO_HEADER* H, bigint n, FILE* output_file);
bigint mda_write_uint16(uint16_t* data, struct MDAIO_HEADER* H, bigint n, FILE* output_file);
bigint mda_write_float64(double* data, struct MDAIO_HEADER* H, bigint n, FILE* output_file);
bigint mda_write_uint32(uint32_t* data, struct MDAIO_HEADER* H, bigint n, FILE* output_file);

//here's an example usage function. See top of file for more info.
void transpose_array(char* infile_path, char* outfile_path);

#endif // MDAIO_H
