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
#ifndef ndarray_h
#define ndarray_h

#include <pybind11/numpy.h>
namespace py = pybind11;

template <typename T>
struct NDArray {
    NDArray(py::array_t<T>& X)
    {
        auto buf = X.request();
        ndim = buf.ndim;
        size = buf.size;
        for (int d = 0; d < ndim; d++) {
            shape.push_back(buf.shape[d]);
        }
        ptr = (T*)buf.ptr;
    }

    long int ndim = 0; // The number of dimensions
    long int size = 0; // The total size of the array (product of dimensions)
    std::vector<int> shape; // The sizes of the dimensions
    T* ptr = 0; // Pointer to the actual data
};

#endif