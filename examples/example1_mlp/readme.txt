mlp-run synthesize_v1.mlp synthesize --samplerate=30000 --duration=600 --timeseries=raw.mda --geom=geom.csv --waveforms_true=waveforms_true.mda --num_channels=10 --num_units=50
mlp-run ../../pipelines/mountainsort3.mlp sort --raw=raw.mda --geom=geom.csv --firings_out=firings2.mda --_params=params.json --curate=true
mountainview --raw=raw.mda --firings=firings2.mda --geom=geom.csv --samplerate=30000
