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
#include "mdaio.h"
#include "usagetracking.h"
#include <vector>
#include <cstring>
#include <inttypes.h>
#include <QDebug>

//can be replaced by std::is_same when C++11 is enabled
template <class T, class U>
struct is_same {
    enum {
        value = 0
    };
};
template <class T>
struct is_same<T, T> {
    enum {
        value = 1
    };
};

int mda_get_num_bytes_per_entry(int data_type)
{
    int num_bytes_per_entry = 0;
    //make sure we have the right number of bytes per entry in case the programmer forgot to set this!
    if (data_type == MDAIO_TYPE_BYTE)
        num_bytes_per_entry = 1;
    else if (data_type == MDAIO_TYPE_COMPLEX)
        num_bytes_per_entry = 8;
    else if (data_type == MDAIO_TYPE_FLOAT32)
        num_bytes_per_entry = 4;
    else if (data_type == MDAIO_TYPE_INT16)
        num_bytes_per_entry = 2;
    else if (data_type == MDAIO_TYPE_INT32)
        num_bytes_per_entry = 4;
    else if (data_type == MDAIO_TYPE_UINT16)
        num_bytes_per_entry = 2;
    else if (data_type == MDAIO_TYPE_FLOAT64)
        num_bytes_per_entry = 8;
    else if (data_type == MDAIO_TYPE_UINT32)
        num_bytes_per_entry = 4;
    return num_bytes_per_entry;
}

bigint mda_read_header(struct MDAIO_HEADER* HH, FILE* input_file)
{
    bigint num_read = 0;
    bigint i;
    size_t totsize;
    bool uses64bitdims = false;

    //initialize
    HH->data_type = 0;
    HH->num_bytes_per_entry = 0;
    HH->num_dims = 0;
    for (i = 0; i < MDAIO_MAX_DIMS; i++)
        HH->dims[i] = 1;
    HH->header_size = 0;

    if (!input_file)
        return 0;

    //data type
    num_read = jfread(&HH->data_type, 4, 1, input_file);
    if (num_read < 1) {
        printf("mda_read_header: Problem reading input file **.\n");
        return 0;
    }

    if ((HH->data_type < -7) || (HH->data_type >= 0)) {
        printf("mda_read_header: Problem with data type:  %d.\n", HH->data_type);
        return 0;
    }

    //number of bytes per entry
    num_read = jfread(&HH->num_bytes_per_entry, 4, 1, input_file);
    if (num_read < 1)
        return 0;

    if ((HH->num_bytes_per_entry <= 0) || (HH->num_bytes_per_entry > 8))
        return 0;

    if (HH->num_bytes_per_entry != mda_get_num_bytes_per_entry(HH->data_type)) {
        qWarning() << "In mda_read_header: num_bytes_per_entry does not agree with data_type. Setting proper value." << HH->num_bytes_per_entry << HH->data_type << mda_get_num_bytes_per_entry(HH->data_type);
        HH->num_bytes_per_entry = mda_get_num_bytes_per_entry(HH->data_type);
    }

    //number of dimensions
    num_read = jfread(&HH->num_dims, 4, 1, input_file);
    if (num_read < 1)
        return 0;

    if ((HH->num_dims == 0) || (HH->num_dims > MDAIO_MAX_DIMS))
        return 0;

    if (HH->num_dims < 0) {
        uses64bitdims = true;
        HH->num_dims = -HH->num_dims;
    }

    //the dimensions
    totsize = 1;
    if (uses64bitdims) {
        for (i = 0; i < HH->num_dims; ++i) {
            uint64_t dim0;
            num_read = jfread(&dim0, sizeof(dim0), 1, input_file);
            if (num_read < 1)
                return 0;
            HH->dims[i] = dim0;
            totsize *= HH->dims[i];
        }
    }
    else {
        for (i = 0; i < HH->num_dims; i++) {
            int32_t dim0;
            num_read = jfread(&dim0, sizeof(dim0), 1, input_file);
            if (num_read < 1)
                return 0;
            if (dim0 < 0) {
                printf("mda_read_header: Dimension %ld less than 0: %d\n", i + 1, dim0);
            }
            HH->dims[i] = dim0;
            totsize *= HH->dims[i];
        }
    }

    if (totsize > MDAIO_MAX_SIZE) {
        printf("mda_read_header: Problem with total size: %ld\n", totsize);
        return 0;
    }

    HH->header_size = 3 * sizeof(int32_t) + HH->num_dims * (uses64bitdims ? sizeof(uint64_t) : sizeof(int32_t));

    //we're done!
    return 1;
}

bigint mda_write_header(struct MDAIO_HEADER* X, FILE* output_file)
{
    bigint num_bytes;
    bigint i;

    if ((X->num_dims <= 0) || (X->num_dims > MDAIO_MAX_DIMS)) {
        printf("mda_write_header: Problem with num dims: %d\n", X->num_dims);
        return 0;
    }
    bool uses64bitdims = false;
    for (int i = 0; !uses64bitdims && (i < X->num_dims); ++i) {
        if (X->dims[i] > 2e9)
            uses64bitdims = true;
    }

    //data type
    num_bytes = fwrite(&X->data_type, 4, 1, output_file);
    if (num_bytes < 1)
        return 0;

    //number of bytes per entry
    X->num_bytes_per_entry = mda_get_num_bytes_per_entry(X->data_type);
    num_bytes = fwrite(&X->num_bytes_per_entry, 4, 1, output_file);
    if (num_bytes < 1)
        return 0;

    //number of dimensions
    int32_t num_dims = uses64bitdims ? -X->num_dims : X->num_dims;
    num_bytes = fwrite(&num_dims, sizeof(num_dims), 1, output_file);
    if (num_bytes < 1)
        return 0;

    //the dimensions
    for (i = 0; i < X->num_dims; i++) {
        if (uses64bitdims) {
            uint64_t dim0 = X->dims[i];
            num_bytes = fwrite(&dim0, sizeof(dim0), 1, output_file);
        }
        else {
            int32_t dim0 = X->dims[i];
            num_bytes = fwrite(&dim0, sizeof(dim0), 1, output_file);
        }
        if (num_bytes < 1)
            return 0;
    }

    X->header_size = 3 * sizeof(int32_t) + X->num_dims * (uses64bitdims ? sizeof(uint64_t) : sizeof(int32_t));

    //we're done!
    return 1;
}

template <typename SourceType, typename TargetType>
bigint mdaReadData_impl(TargetType* data, const bigint size, FILE* inputFile)
{
    if (is_same<TargetType, SourceType>::value) {
        return jfread(data, sizeof(SourceType), size, inputFile);
    }
    else {
        std::vector<SourceType> tmp(size);
        const bigint ret = jfread(&tmp[0], sizeof(SourceType), size, inputFile);
        std::copy(tmp.begin(), tmp.end(), data);
        return ret;
    }
}

template <typename Type>
bigint mdaReadData(Type* data, const struct MDAIO_HEADER* header, const bigint size, FILE* inputFile)
{
    if (header->data_type == MDAIO_TYPE_BYTE) {
        return mdaReadData_impl<unsigned char>(data, size, inputFile);
    }
    else if (header->data_type == MDAIO_TYPE_FLOAT32) {
        return mdaReadData_impl<float>(data, size, inputFile);
    }
    else if (header->data_type == MDAIO_TYPE_INT16) {
        return mdaReadData_impl<int16_t>(data, size, inputFile);
    }
    else if (header->data_type == MDAIO_TYPE_INT32) {
        return mdaReadData_impl<int32_t>(data, size, inputFile);
    }
    else if (header->data_type == MDAIO_TYPE_UINT16) {
        return mdaReadData_impl<uint16_t>(data, size, inputFile);
    }
    else if (header->data_type == MDAIO_TYPE_FLOAT64) {
        return mdaReadData_impl<double>(data, size, inputFile);
    }
    else if (header->data_type == MDAIO_TYPE_UINT32) {
        return mdaReadData_impl<uint32_t>(data, size, inputFile);
    }
    else
        return 0;
}

template <typename TargetType, typename DataType>
bigint mdaWriteData_impl(DataType* data, const bigint size, FILE* outputFile)
{
    if (is_same<DataType, TargetType>::value) {
        return fwrite(data, sizeof(DataType), size, outputFile);
    }
    else {
        std::vector<TargetType> tmp(size);
        std::copy(data, data + size, tmp.begin());
        return fwrite(&tmp[0], sizeof(TargetType), size, outputFile);
    }
}

template <typename DataType>
bigint mdaWriteData(DataType* data, const bigint size, const struct MDAIO_HEADER* header, FILE* outputFile)
{
    if (header->data_type == MDAIO_TYPE_BYTE) {
        return mdaWriteData_impl<unsigned char>(data, size, outputFile);
    }
    else if (header->data_type == MDAIO_TYPE_FLOAT32) {
        return mdaWriteData_impl<float>(data, size, outputFile);
    }
    else if (header->data_type == MDAIO_TYPE_INT16) {
        return mdaWriteData_impl<int16_t>(data, size, outputFile);
    }
    else if (header->data_type == MDAIO_TYPE_INT32) {
        return mdaWriteData_impl<int32_t>(data, size, outputFile);
    }
    else if (header->data_type == MDAIO_TYPE_UINT16) {
        return mdaWriteData_impl<uint16_t>(data, size, outputFile);
    }
    else if (header->data_type == MDAIO_TYPE_FLOAT64) {
        return mdaWriteData_impl<double>(data, size, outputFile);
    }
    else if (header->data_type == MDAIO_TYPE_UINT32) {
        return mdaWriteData_impl<uint32_t>(data, size, outputFile);
    }
    else
        return 0;
}

bigint mda_read_byte(unsigned char* data, struct MDAIO_HEADER* H, bigint n, FILE* input_file)
{
    return mdaReadData(data, H, n, input_file);
}

bigint mda_read_float32(float* data, struct MDAIO_HEADER* H, bigint n, FILE* input_file)
{
    return mdaReadData(data, H, n, input_file);
}

bigint mda_read_float64(double* data, struct MDAIO_HEADER* H, bigint n, FILE* input_file)
{
    return mdaReadData(data, H, n, input_file);
}

bigint mda_read_int16(int16_t* data, struct MDAIO_HEADER* H, bigint n, FILE* input_file)
{
    return mdaReadData(data, H, n, input_file);
}

bigint mda_read_int32(int32_t* data, struct MDAIO_HEADER* H, bigint n, FILE* input_file)
{
    return mdaReadData(data, H, n, input_file);
}

bigint mda_read_uint16(uint16_t* data, struct MDAIO_HEADER* H, bigint n, FILE* input_file)
{
    return mdaReadData(data, H, n, input_file);
}

bigint mda_read_uint32(uint32_t* data, struct MDAIO_HEADER* H, bigint n, FILE* input_file)
{
    return mdaReadData(data, H, n, input_file);
}

bigint mda_write_byte(unsigned char* data, struct MDAIO_HEADER* H, bigint n, FILE* output_file)
{
    return mdaWriteData(data, n, H, output_file);
}

bigint mda_write_float32(const float* data, struct MDAIO_HEADER* H, bigint n, FILE* output_file)
{
    return mdaWriteData(data, n, H, output_file);
}

bigint mda_write_int16(int16_t* data, struct MDAIO_HEADER* H, bigint n, FILE* output_file)
{
    return mdaWriteData(data, n, H, output_file);
}

bigint mda_write_int32(int32_t* data, struct MDAIO_HEADER* H, bigint n, FILE* output_file)
{
    return mdaWriteData(data, n, H, output_file);
}

bigint mda_write_uint16(uint16_t* data, struct MDAIO_HEADER* H, bigint n, FILE* output_file)
{
    return mdaWriteData(data, n, H, output_file);
}

bigint mda_write_float64(double* data, struct MDAIO_HEADER* H, bigint n, FILE* output_file)
{
    return mdaWriteData(data, n, H, output_file);
}

bigint mda_write_uint32(uint32_t* data, struct MDAIO_HEADER* H, bigint n, FILE* output_file)
{
    return mdaWriteData(data, n, H, output_file);
}

void mda_copy_header(struct MDAIO_HEADER* ret, const struct MDAIO_HEADER* X)
{
    std::memcpy(ret, X, sizeof(*ret));
}

void transpose_array(char* infile_path, char* outfile_path)
{
    //open the files for reading/writing and declare some variables
    FILE* infile = fopen(infile_path, "rb");
    FILE* outfile = fopen(outfile_path, "wb");
    struct MDAIO_HEADER H;
    bigint M, N;
    bigint i, j;
    float* data_in, *data_out;

    if (!infile)
        return;
    if (!outfile) {
        fclose(infile);
        return;
    }

    //read the header
    mda_read_header(&H, infile);

    //if the data type is zero then there was a problem reading
    if (!H.data_type) {
        fclose(infile);
        fclose(outfile);
        return;
    }

    //get the dimensions and allocate the in/out arrays
    M = H.dims[0];
    N = H.dims[1];
    data_in = (float*)malloc(sizeof(float) * M * N);
    data_out = (float*)malloc(sizeof(float) * M * N);

    //Read the data -- note that we don't care what the actual type is.
    //This is a great trick!
    //See top of file for more info
    //Note that we could have decided to read only portions of the file if
    //N*M is too large for memory
    mda_read_float32(data_in, &H, M * N, infile);

    //Perform the transpose
    for (j = 0; j < N; j++) {
        for (i = 0; i < M; i++) {
            data_out[i + N * j] = data_in[j + M * i];
        }
    }

    //Swap the dimensions and write the output data
    H.dims[0] = N;
    H.dims[1] = M;
    mda_write_float32(data_out, &H, M * N, outfile);

    //clean up and we're done
    free(data_in);
    free(data_out);
    fclose(infile);
    fclose(outfile);
}
