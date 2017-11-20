FROM ubuntu:xenial

MAINTAINER Jeremy Magland

RUN apt-get update

# Install Qt5
RUN apt-get update && apt-get install -y software-properties-common
RUN apt-add-repository -y ppa:ubuntu-sdk-team/ppa
RUN apt-get update && apt-get install -y qtdeclarative5-dev
RUN apt-get update && apt-get install -y qt5-default
RUN apt-get update && apt-get install -y qtbase5-dev
RUN apt-get update && apt-get install -y qtscript5-dev

RUN apt-get update && apt-get install -y libqt5webkit5-dev

# Install
RUN apt-get update && apt-get install -y nodejs npm nodejs-legacy

RUN apt-get update && apt-get install -y python3 python3-pip
RUN pip3 install numpy scipy pybind11 cppimport numpydoc

RUN apt-get update && apt-get install -y libfftw3-dev
RUN pip3 install sklearn

RUN apt-get update && apt-get install -y git

RUN mkdir /work && echo "4"
WORKDIR /work
RUN git clone https://github.com/flatironinstitute/mountainlab
RUN git clone https://github.com/flatironinstitute/mlpipeline
RUN git clone https://github.com/flatironinstitute/mountainview

WORKDIR /work/mountainlab 
RUN ./compile_components.sh
ENV PATH="/work/mountainlab/bin:${PATH}"

WORKDIR /work/mlpipeline
RUN ./compile_components.sh
ENV PATH="/work/mlpipeline/bin:${PATH}"
ENV PATH="/work/mlpipeline/utils/mlp:${PATH}"

WORKDIR /work/mountainview
RUN ./compile_components.sh
ENV PATH="/work/mountainview/bin:${PATH}"


#### MountainSort 
WORKDIR /work/mountainlab/packages
RUN git clone https://github.com/flatironinstitute/mountainsort && echo "3"

WORKDIR /work/mountainlab/packages/mountainsort
RUN ./compile_components.sh

RUN mp-list-processors

#### Run example
WORKDIR /work
RUN git clone https://github.com/flatironinstitute/mountainsort_examples && echo "3"

WORKDIR /work/mountainsort_examples/devel/accuracy_test_1

RUN ./run_test.sh

