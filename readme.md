ciglet
===

Ciglet is a lightweight C library for digital signal processing, in particular audio and speech processing.

Language: C99

License: BSD

Background
---

In the past few years I've been writing and maintaining quite a few C-written speech processing projects (e.g. *libllsm*, *libpyin*, *moresampler*), some of which were partially translated from Matlab/Octave-based prototypes. Gradually I ended up having lots of frequent rewrites of C version of Matlab routines such as `sum`, `fir1`, `conv`, `interp1`, etc. In October 2016 I finally made the move to extract all these repeating patterns and group them into one independent library. Since it would be referred by quite a few projects, the new library is meant to be lightweight, easy to link, and fast, which is why it's named *ciglet*, an acronym of "C-written sIGnal codeLETs".

Function index
---

### Scalar operations

`max`, `min`, `linterp`, `fastatan2`, `randu`, `randn`

### Complex scalar operations

`c_cplx`, `c_conj`, `c_add`, `c_sub`, `c_mul`, `c_div`, `c_exp`, `c_abs`, `c_arg`

### Vector operations and statistics

`sumfp`, `sumsqrfp`, `maxfp`, `minfp`, `meanfp`, `varfp`, `medianfp`, `selectnth`, `sort`, `xcorr`, `corr`, `cov`, `find_peak`/`find_valley`

### Numerical routines

`fzero`, `polyval`, `roots`

### Memory (de)allocation

`linspace`, `iota`, `malloc2d`, `free2d`, `copy2d`, `flatten`, `reshape`, `transpose`

### Audio I/O

`wavread`, `wavwrite`

### General DSP routines

`fetch_frame`, `gensin`/`gensins`, `boxcar`, `hanning`, `hamming`, `mltsine`, `blackman_harris`, `blackman`, `fft`/`ifft`, `idft`, `dct`, `fftshift`, `wrap`/`unwrap`, `diff`, `cumsum`, `flip`, `abscplx`, `argcplx`, `polar2real`, `polar2imag`, `phase_diff`, `complete_symm`/`complete_asymm`, `rceps`, `irceps`, `minphase`, `fir1`, `conv`, `filter`, `filtfilt`, `levinson`, `lpc`/`flpc`, `lpgain`, `lpspec`, `lpresf`, `interp1`/`interp1u`, `sincinterp1u`, `medfilt1`, `white_noise`, `moving_avg`, `moving_rms`, `itakura_saito`, `safe_aliased_sinc`/`safe_aliased_dsinc`, `rresample`

### Audio/speech processing routines

`mel2freq`/`freq2mel`, `freq2bark`/`bark2freq`, `eqloud`, `melspace`, `correlogram`, `invcrgm`, `stft`/`istft`, `qifft`, `spgm2cegm`/`cegm2spgm`, `create_filterbank`, `delete_filterbank`, `filterbank_spgm`, `filterbank_spec`, `be2cc`, `be2ccgm`, `spec2env`, `lfmodel_from_rd`, `lfmodel_spectrum`, `lfmodel_period`

### Plotting utilities (Gnuplot interface, unavailable on Windows)

`plotopen`, `plot`, `imagesc`, `plotclose`
