FROM ubuntu:16.04

MAINTAINER Jeremy Magland

# Python3 with numpy
RUN apt-get update && \
    apt-get install -y \
    python3 python3-numpy python3-pip

RUN pip3 install pybind11 cppimport
RUN pip3 install scipy
RUN pip3 install sklearn

RUN apt-get update && apt-get install -y fftw3-dev

ADD . /package

# Build
WORKDIR /package
RUN basic/basic.mp spec > basic.spec
RUN drift/drift.mp spec > drift.spec
RUN synthesis/synthesis.mp spec > synthesis.spec
