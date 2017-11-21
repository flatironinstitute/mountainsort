# Set the base image to Ubuntu
FROM ubuntu:xenial

MAINTAINER Jeremy Magland

RUN apt-get update && echo "1"

RUN apt-get update && apt-get install -y software-properties-common

RUN add-apt-repository -y ppa:magland/mountainlab && echo "7"
RUN apt-get update && apt-get install -y mountainlab
RUN apt-get update && apt-get install -y mlpipeline
RUN apt-get update && apt-get install -y mountainsort && echo "1"
RUN apt-get update && apt-get install -y mountainview

RUN apt-get update && apt-get install -y python3-pip
RUN pip3 install --upgrade numpy

RUN mp-list-processors

RUN apt-get update && apt-get install -y git

#### Run example
WORKDIR /work
RUN git clone https://github.com/flatironinstitute/mountainsort_examples && echo "2"

# WORKDIR /work/mountainsort_examples/example1_mlp

# RUN mlp-run synthesize_v1.mlp synthesize --samplerate=30000 --duration=600 --timeseries=raw.mda --geom=geom.csv --waveforms_true=waveforms_true.mda --num_channels=10 --num_units=50
# RUN mlp-run mountainsort3.mlp sort --raw=raw.mda --geom=geom.csv --firings_out=firings.mda --_params=params.json

# RUN mda firings.mda

WORKDIR /work/mountainsort_examples/devel/accuracy_test_1
RUN ./run_test.sh



