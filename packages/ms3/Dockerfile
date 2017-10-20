# Set the base image to Ubuntu
FROM ubuntu:16.04

MAINTAINER Jeremy Magland

RUN apt-get update

# Install qt5
RUN apt-get update && apt-get install -y software-properties-common
RUN apt-get update && apt-add-repository ppa:ubuntu-sdk-team/ppa
RUN apt-get update && apt-get install -y qtdeclarative5-dev qt5-default qtbase5-dev qtscript5-dev
RUN apt-get update && apt-get install -y make g++

# Install fftw3
RUN apt-get update && apt-get install -y fftw3-dev

ADD . /package

# Build
WORKDIR /package
RUN qmake
RUN make -j
RUN bin/ms3.mp spec > ms3.spec
