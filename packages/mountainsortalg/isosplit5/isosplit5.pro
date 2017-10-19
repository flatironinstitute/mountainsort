QT += core
QT -= gui

CONFIG += c++11

DESTDIR = bin
OBJECTS_DIR = build
TARGET = isosplit5-dummy
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp \
    isocut5.cpp \
    jisotonic5.cpp \
    isosplit5.cpp

HEADERS += \
    isocut5.h \
    jisotonic5.h \
    isosplit5.h
