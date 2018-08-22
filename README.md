# spectral_analysis
spectral analysis library in matlab

For computing spectrograms:
- **mtgram**: computes a spectrogram with pmtm instead of periodogram
- **pchavegram**: computes a spectrogram with pchave instead of periodogram (unfinished)

For processing time series:
- **prewhiten**: prewhitens a time series by fitting an AR(p) model and inverse filtering the time series to obtain a whitened version of the series.

TO DO:
- [ ] Implement pchavegram (just a copy of mtgram right now)
- [x] Generalize prewhiten for any FIR filters
