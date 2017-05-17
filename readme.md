ciglet
===

Ciglet is a lightweight C library for digital signal processing, in particular audio and speech processing.

Language: C99

License: BSD

Background
---

In the past few years I've been writing and maintaining quite a few C-written speech processing projects (e.g. *libllsm*, *libpyin*, *moresampler*), some of which were partially translated from Matlab/Octave-based prototypes. Gradually I ended up having lots of frequent rewrites of C version of Matlab routines such as `sum`, `fir1`, `conv`, `interp1`, etc. In August 2016 I finally made the move to extract all these repeating patterns and group them into one independent library. Since it would be referred by quite a few projects, the new library is meant to be lightweight, easy to link, and fast, which is why it's named *ciglet*, an acronym of "C-written sIGnal codeLETs".

Function index
---

### Scalar operations

* random number generation: `randu`, `randn`
* miscellaneous: `max`, `min`, `linterp`, `fastatan2`
* complex arithmetics: `c_cplx`, `c_conj`, `c_add`, `c_sub`, `c_mul`, `c_div`, `c_exp`, `c_abs`, `c_arg`

### Vector operations and statistics

* vectorized arithmetics: `sumfp`, `sumsqrfp`, `maxfp`, `minfp`
* descriptive statistics: `meanfp`, `varfp`, `medianfp`, `xcorr`, `corr`, `cov`
* sorting: `selectnth`, `sort`
* peak picking: `find_peak`, `find_valley`, `find_maxima`, `find_minima`

### Numerical routines

`fzero`, `polyval`, `roots`

### Basic linear algebra

* products: `matmul`, `mvecmul`, `dot`
* solving a linear system: `lu`, `lusolve`
* pivoting: `ppivot`, `permm`, `permv`

### Memory (de)allocation

* enumeration: `linspace`, `iota`
* 2d array operations: `malloc2d`, `free2d`, `copy2d`, `flatten`, `reshape`, `transpose`

### Audio I/O

`wavread`, `wavwrite`

### General DSP routines

* windows: `boxcar`, `hanning`, `hamming`, `mltsine`, `blackman_harris`, `nuttall98`, `blackman`
* Fourier transform: `fft`, `ifft`, `czt`, `iczt`, `idft`, `dct`, `fftshift`
* phase manipulation: `wrap`, `unwrap`, `phase_diff`
* complex number conversion: `abscplx`, `argcplx`, `polar2real`, `polar2imag`, `complete_symm`, `complete_asymm`
* cepstral analysis: `rceps`, `irceps`, `minphase`
* filtering: `fir1`, `conv`, `filter`, `filtfilt`, `moving_avg`, `moving_rms`, `medfilt1`, `kalmanf1d`, `kalmans1d`
* linear prediction: `levinson`, `lpc`, `flpc`, `lpgain`, `lpspec`, `lpresf`
* interpolation: `interp1`, `interp1u`, `sincinterp1u`, `interp_in_blank`, `rresample`
* operations on sinusoids: `gensin`, `gensins`, `safe_aliased_sinc`, `safe_aliased_dsinc`
* miscellaneous: `fetch_frame`, `diff`, `cumsum`, `flip` ,`white_noise`, `itakura_saito`

### Audio/speech processing routines

* psychoacoustics: `mel2freq`, `freq2mel`, `freq2bark`, `bark2freq`, `eqloud`, `melspace`
* frequency estimation: `ifdetector_estimate`, `correlogram`, `invcrgm`
* spectrogram and STFT: `stft`, `istft`, `qifft`, `spgm2cegm`, `cegm2spgm`
* filterbank analysis: `filterbank_spgm`, `filterbank_spec`, `be2cc`, `be2ccgm`
* spectral envelope estimation: `spec2env`
* glottal model: `lfmodel_from_rd`, `lfmodel_spectrum`, `lfmodel_period`

### Plotting utilities (Gnuplot interface, unavailable on Windows)

`plotopen`, `plot`, `imagesc`, `plotclose`
