mlp-run synthesize_v1.mlp synthesize --duration=60 --samplerate=30000 --timeseries=raw.mda
mlp-run ../../pipelines/mountainsort3.mlp sort --firings_out=firings.mda --raw=raw.mda
mountainview --raw=raw.mda --firings=firings.mda --samplerate=30000
