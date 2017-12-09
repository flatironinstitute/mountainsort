def WMLDR(data, wname="db4", maxlevel=6, mode='sym'):
    """ Function by Martin Spacek from https://github.com/spyke/spyke

    Perform wavelet multi-level decomposition and reconstruction (WMLDR) on multichannel
    data. See Wiltschko2008. Default to Daubechies(4) wavelet. Modifies data in-place, at
    least for now. The effective cutoff frequency is:

    fc = (sampfreq / 2) / 2**maxlevel                     (Wiltschko2008)

    For sampfreq of 25 kHz and maxlevel of 6, the effective cutoff frequency is 195 Hz.
    For sampfreq of 30 kHz and maxlevel of 6, the effective cutoff frequency is 234 Hz.

    TODO: for now, this only returns highpass data. In the future, this probably should
    return both low and highpass data (and not modify it in-place). The Discussion in
    Wiltschko2008 suggests that this approach cannot be used to extract the LFP, but
    I don't see why you can't simply subtract the highpass data from the raw data to get the
    lowpass data.

    Signal extension modes (from PyWavelets docs):

    PyWavelets provides several methods of signal extrapolation that can be used to minimize
    edge effects. PyWavelet's default is 'sym':

    zpd - zero-padding - signal is extended by adding zero samples:
    ... 0  0 | x1 x2 ... xn | 0  0 ...

    cpd - constant-padding - border values are replicated:
    ... x1 x1 | x1 x2 ... xn | xn xn ...

    sym - symmetric-padding - signal is extended by mirroring samples:
    ... x2 x1 | x1 x2 ... xn | xn xn-1 ...

    ppd - periodic-padding - signal is treated as a periodic one:
    ... xn-1 xn | x1 x2 ... xn | x1 x2 ...

    sp1 - smooth-padding - signal is extended according to the first derivatives calculated on
    the edges (straight line)

    DWT performed for these extension modes is slightly redundant, but ensures perfect
    reconstruction. To receive the smallest possible number of coefficients, computations can
    be performed with the periodization mode:

    per - periodization - is like periodic-padding but gives the smallest possible number of
    decomposition coefficients. IDWT must be performed with the same mode.
    """
    import pywt

    data = np.atleast_2d(data)
    nt = data.shape[1]
    # reconstructed signals always seem to have an even number of points. If the number of
    # input data points is odd, trim last data point from reconstructed signal:
    isodd = nt % 2
    # filter data in place, iterate over channels in rows:
    nchans = len(data)
    for chani in range(nchans):
        # decompose the signal:
        cs = pywt.wavedec(data[chani], wname, mode=mode, level=maxlevel)
        # destroy the appropriate approximation coefficients to get highpass data:
        cs[0] = None
        # reconstruct the signal:
        recsignal = pywt.waverec(cs, wname, mode=mode)
        ntrec = len(recsignal)
        data[chani] = recsignal[:ntrec-isodd]

    return data
