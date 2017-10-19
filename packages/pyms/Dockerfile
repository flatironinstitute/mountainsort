# Set the base image to Ubuntu
FROM ubuntu:16.04

MAINTAINER Jeremy Magland

RUN apt-get update && apt-get install -y python3 python3-pip

RUN pip3 install numpy pybind11 cppimport scipy

RUN apt-get update && apt-get install -y fftw3-dev

ADD . /package

# Build
WORKDIR /package
RUN basic/basic.mp spec > basic.spec
RUN drift/drift.mp spec > drift.spec
RUN synthesis/synthesis.mp spec > synthesis.spec
