#!/bin/bash

set -e

duration=600
samplerate=30000
waveform_upsamplefac=13

# Synthesize some raw data (create waveforms.mda, geom.csv and raw.mda)
mp-run-process pyms.synthesize_random_waveforms \
			--waveforms_out=waveforms.mda --geometry_out=geom.csv --upsamplefac=$waveform_upsamplefac

mp-run-process pyms.synthesize_random_firings \
			--firings_out=firings_true.mda --samplerate=$samplerate --duration=$duration

mp-run-process pyms.synthesize_timeseries \
			--waveforms=waveforms.mda --firings=firings_true.mda \
			--timeseries_out=raw.mda \
			--waveform_upsamplefac=$waveform_upsamplefac \
			--samplerate=$samplerate --duration=$duration

# Preprocessing (create filt.mda and pre.mda)
mp-run-process pyms.bandpass_filter \
			--timeseries=raw.mda --timeseries_out=filt.mda \
			--samplerate=$samplerate --freq_min=300 --freq_max=6000 --freq_wid=1000

mp-run-process pyms.normalize_channels \
			--timeseries=filt.mda --timeseries_out=pre.mda

# Sorting (create firings.mda)
mp-run-process mountainsortalg.ms3 \
			--timeseries=pre.mda --firings_out=firings.mda \
			--detect_threshold=3 --detect_sign=1

# View the results (launch the viewer)
mountainview --raw=raw.mda --filt=filt.mda --pre=pre.mda --firings=firings.mda --samplerate=$samplerate