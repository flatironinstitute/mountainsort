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

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "ndarray.h"


namespace py = pybind11;

// A custom error handler to propagate errors back to Python
class error : public std::exception {
public:
    error(const std::string& msg)
        : msg_(msg){};
    virtual const char* what() const throw()
    {
        return msg_.c_str();
    }

private:
    std::string msg_;
};

void extract_clips(int M,int T,int L,int N,float *clips_out,float *X,int *times,int clip_size) {
    int Tmid=(int)((T+1)/2-1);
    for (int i=0; i<L; i++) {
        int t=times[i];
        if ((0<=t-Tmid)&&(t-Tmid+T-1<N)) {
            int ii=M*(t-Tmid);
            int jj=M*T*i;
            for (int m=0; m<M; m++) {
                for (int a=0; a<T; a++) {
                    clips_out[jj]=X[ii];
                    ii++;
                    jj++;
                }
            }
        }
    }
}

void extract_clips_(py::array_t<float> &clips_out,py::array_t<float> X, py::array_t<int> times, int clip_size)
{
    NDArray<float> Xa(X);
    NDArray<float> Ca(clips_out);
    NDArray<int> Ta(times);
    if (Xa.ndim != 2) {
        printf ("Incorrect number of dimensions: %ld\n", Xa.ndim);
        throw error("Incorrect number of dimensions");
    }
    int L = Ta.size;
    int M = Xa.shape[0];
    int N = Xa.shape[1];
    int T = clip_size;

    if ((Ca.shape[0]!=M)||(Ca.shape[1]!=T)||(Ca.shape[2]!=L)) {
        printf ("%d,%d,%d  <>  %d,%d,%d",Ca.shape[0],Ca.shape[1],Ca.shape[2],M,T,L);
        throw error("Incorrect size of clips_out");
    }

    extract_clips(M,T,L,N,Ca.ptr,Xa.ptr,Ta.ptr,clip_size);
}



PYBIND11_MODULE(basic_cpp,m)
{
    m.def("extract_clips", &extract_clips_, "", py::arg().noconvert(),py::arg(),py::arg(),py::arg());
}

/*
PYBIND11_PLUGIN(basic_cpp) {
    py::module m("basic_cpp");
    m.def("extract_clips", &extract_clips_, "", py::arg().noconvert(),py::arg(),py::arg(),py::arg());
    return m.ptr();
}
*/
