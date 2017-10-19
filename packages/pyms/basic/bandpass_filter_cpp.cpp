/*
<%
cfg['include_dirs'] = ['../mlpy']
cfg['sources'] = ['bandpass_filter_kernel.cpp']
cfg['linker_args'] = ['-fopenmp','-lfftw3']
cfg['compiler_args'] = ['-fopenmp']
setup_pybind11(cfg)
%>



*/
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "ndarray.h"
#include "bandpass_filter_kernel.h"

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

void bandpass_filter_(py::array_t<float> &X,double samplerate,double freq_min,double freq_max,double freq_wid)
{
    NDArray<float> Xa(X);
    if (Xa.ndim != 2) {
        printf ("Incorrect number of dimensions: %ld\n", Xa.ndim);
        throw error("Incorrect number of dimensions");
    }
    int M = Xa.shape[0];
    int N = Xa.shape[1];

    //bandpass_filter_kernel(M,N,Xa.ptr,samplerate,freq_min,freq_max,freq_wid,0,N-1);
    bandpass_filter_kernel_multithread(M,N,Xa.ptr,samplerate,freq_min,freq_max,freq_wid);
}

PYBIND11_MODULE(bandpass_filter_cpp,m)
{
    m.def("bandpass_filter", &bandpass_filter_, "", py::arg().noconvert(),py::arg(),py::arg(),py::arg(),py::arg());
}


