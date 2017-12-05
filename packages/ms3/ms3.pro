QT += core
QT -= gui
QT += qml

CONFIG += c++11

DESTDIR = bin
OBJECTS_DIR = build
MOC_DIR= build
TARGET = ms3.mp
TEMPLATE = app

QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS += -fopenmp

#FFTW
USE_FFTW3=$$(USE_FFTW3)
CONFIG("no_fftw3") {
    warning(Not using FFTW3)
}
else {
    LIBS += -lfftw3f -lfftw3f_threads
    SOURCES += p_spikeview_templates.cpp
    HEADERS += p_spikeview_templates.h
}

#-std=c++11   # AHB removed since not in GNU gcc 4.6.3

SOURCES += \
    mountainsort3_main.cpp \
    dimension_reduce_clips.h \
    dimension_reduce_clips.cpp \
    reorder_labels.cpp \
    p_run_metrics_script.cpp \
    p_spikeview_metrics.cpp \
    p_synthesize_timeseries.cpp \
    p_combine_firing_segments.cpp \
    p_extract_firings.cpp \
    p_concat_timeseries.cpp \
    p_banjoview_cross_correlograms.cpp \
    kdtree.cpp \
    get_sort_indices.cpp \
    p_create_multiscale_timeseries.cpp \
    p_mv_discrimhist.cpp

HEADERS += \
    mountainsort3_main.h \
    reorder_labels.h \
    p_run_metrics_script.h \
    p_spikeview_metrics.h \
    p_synthesize_timeseries.h \
    p_combine_firing_segments.h \
    p_extract_firings.h \
    kdtree.h \
    get_sort_indices.h \
    p_create_multiscale_timeseries.h \
    p_mv_discrimhist.h

HEADERS += \
    p_extract_clips.h \
    p_create_firings.h \
    p_combine_firings.h \
    p_whiten.h \
    p_apply_timestamp_offset.h \
    p_link_segments.h \
    p_cluster_metrics.h \
    p_split_firings.h \
    p_concat_firings.h \
    p_compute_templates.h \
    p_load_test.h \
    p_compute_amplitudes.h \
    p_isolation_metrics.h \
    p_confusion_matrix.h \
    p_reorder_labels.h \
    p_bandpass_filter.h \
    p_mask_out_artifacts.h \
    p_mv_compute_templates.h \
    p_mv_compute_amplitudes.h extract_clips.h get_principal_components.h \
    hungarian.h

SOURCES += \
    p_extract_clips.cpp \
    p_create_firings.cpp \
    p_combine_firings.cpp \
    p_whiten.cpp \
    p_apply_timestamp_offset.cpp \
    p_link_segments.cpp \
    p_cluster_metrics.cpp \
    p_split_firings.cpp \
    p_concat_firings.cpp \
    p_compute_templates.cpp \
    p_load_test.cpp \
    p_compute_amplitudes.cpp \
    p_isolation_metrics.cpp \
    p_confusion_matrix.cpp \
    p_reorder_labels.cpp \
    p_bandpass_filter.cpp \
    p_mask_out_artifacts.cpp \
    p_mv_compute_templates.cpp \
    p_mv_compute_amplitudes.cpp extract_clips.cpp get_principal_components.cpp \
    hungarian.cpp

HEADERS += pca.h compute_templates_0.h
SOURCES += pca.cpp compute_templates_0.cpp

INCLUDEPATH += mda
HEADERS += mda/mda_p.h \
    mda/mda.h \
    mda/mda32.h \
    mda/mdaio.h \
    mda/mdareader_p.h \
    mda/mdareader.h \
    mda/usagetracking.h

SOURCES += mda/diskreadmda.cpp \
    mda/diskreadmda32.cpp \
    mda/diskwritemda.cpp \
    mda/mda.cpp \
    mda/mda32.cpp \
    mda/mdaio.cpp \
    mda/mdareader.cpp \
    mda/usagetracking.cpp

INCLUDEPATH += common
HEADERS += common/mlcompute.h \
    common/clparams.h \
    common/textfile.h \
    common/mlutil.h

SOURCES += common/mlcompute.cpp \
    common/clparams.cpp \
    common/textfile.cpp \
    common/mlutil.cpp


