QT += core
QT -= gui

CONFIG += c++11

DESTDIR = bin
OBJECTS_DIR = build
MOC_DIR=build
TARGET = mountainsortalg.mp
TEMPLATE = app

DISTFILES +=

INCLUDEPATH += mda
INCLUDEPATH += isosplit5
INCLUDEPATH += common

QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS += -fopenmp

HEADERS += \
    compute_templates_0.h \
    consolidate_clusters.h \
    detect_events.h \
    fit_stage.h \
    globaltemplatecomputer.h \
    merge_across_channels.h \
    neighborhoodsorter.h \
    p_mountainsort3.h \
    pca.h \
    sort_clips.h \
    mda/diskreadmda.h \
    mda/diskreadmda32.h \
    mda/diskwritemda.h \
    mda/mda_p.h \
    mda/mda.h \
    mda/mda32.h \
    mda/mdaio.h \
    mda/mdareader_p.h \
    mda/mdareader.h \
    mda/usagetracking.h \
    get_sort_indices.h \
    isosplit5/isocut5.h \
    isosplit5/isosplit5.h \
    isosplit5/jisotonic5.h \
    fitstagecomputer.h \
    mountainsortalg_main.h \
    common/mlcompute.h \
    common/clparams.h \
    common/textfile.h \
    common/mlutil.h

SOURCES += \
    compute_templates_0.cpp \
    consolidate_clusters.cpp \
    detect_events.cpp \
    fit_stage.cpp \
    globaltemplatecomputer.cpp \
    merge_across_channels.cpp \
    neighborhoodsorter.cpp \
    p_mountainsort3.cpp \
    pca.cpp \
    sort_clips.cpp \
    mda/diskreadmda.cpp \
    mda/diskreadmda32.cpp \
    mda/diskwritemda.cpp \
    mda/mda.cpp \
    mda/mda32.cpp \
    mda/mdaio.cpp \
    mda/mdareader.cpp \
    mda/usagetracking.cpp \
    get_sort_indices.cpp \
    isosplit5/isocut5.cpp \
    isosplit5/isosplit5.cpp \
    isosplit5/jisotonic5.cpp \
    mountainsortalg_main.cpp \
    fitstagecomputer.cpp \
    common/mlcompute.cpp \
    common/clparams.cpp \
    common/textfile.cpp \
    common/mlutil.cpp

