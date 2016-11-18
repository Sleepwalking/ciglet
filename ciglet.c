#include "ciglet.h"

// fftsg
void cdft(int n, int isgn, FP_TYPE* a);
void rdft(int n, int isgn, FP_TYPE* a);

// fast median
FP_TYPE opt_med3(FP_TYPE* x);
FP_TYPE opt_med5(FP_TYPE* x);
FP_TYPE opt_med7(FP_TYPE* x);
FP_TYPE opt_med9(FP_TYPE* x);
FP_TYPE opt_med25(FP_TYPE* x);

#define def_transpose(nbit, type) \
  static type** cig_transpose_##nbit(type** ptr, size_t m, size_t n) { \
    type** ret = malloc(n * sizeof(type*)); \
    for(size_t i = 0; i < n; i ++) { \
      ret[i] = malloc(m * sizeof(type)); \
      for(size_t j = 0; j < m; j ++) \
        ret[i][j] = ptr[j][i]; \
    } \
    return ret; \
  }

def_transpose(8, unsigned char);
def_transpose(16, uint16_t);
def_transpose(32, uint32_t);
def_transpose(64, uint64_t);

void** cig_transpose(void** ptr, size_t m, size_t n, size_t size) {
  switch(size * 8) {
    case 8:
      return (void**)cig_transpose_8((unsigned char**)ptr, m, n);
      break;
    case 16:
      return (void**)cig_transpose_16((uint16_t**)ptr, m, n);
      break;
    case 32:
      return (void**)cig_transpose_32((uint32_t**)ptr, m, n);
      break;
    case 64:
      return (void**)cig_transpose_64((uint64_t**)ptr, m, n);
      break;
  }
  fprintf(stderr, "Error: invalid element size '%d'.\n", (int)size);
  return NULL;
}

// https://wikicoding.org/wiki/c/Quickselect/
FP_TYPE cig_qselect(FP_TYPE* v, int len, int k) {
# define swap(a, b) { \
    tmp = v[a]; v[a] = v[b]; v[b] = tmp; \
  }
  int i, st;
  FP_TYPE tmp;

  for(st = i = 0; i < len - 1; i ++) {
    if(v[i] > v[len - 1]) continue;
    swap(i, st);
    st ++;
  }
  swap(len - 1, st);

  return k == st ? v[st]
    :st > k ? cig_qselect(v, st, k)
      : cig_qselect(v + st, len - st, k - st);
# undef swap
}

typedef struct {
  FP_TYPE val;
  int idx;
} sort_struct;

int compare_sort_struct(const void* a, const void* b) {
 return ((sort_struct*)a) -> val > ((sort_struct*)b) -> val ? 1 : -1;
}

FP_TYPE* cig_sort(FP_TYPE* x, int nx, int* outidx) {
  sort_struct* s = malloc(nx * sizeof(sort_struct));
  for(int i = 0; i < nx; i ++) {
    s[i].val = x[i];
    s[i].idx = i;
  }
  qsort(s, nx, sizeof(sort_struct), compare_sort_struct);
  FP_TYPE* y = malloc(nx * sizeof(FP_TYPE));
  for(int i = 0; i < nx; i ++)
    y[i] = s[i].val;
  if(outidx != NULL)
    for(int i = 0; i < nx; i ++)
      outidx[i] = s[i].idx;
  free(s);
  return y;
}

FP_TYPE cig_medianfp(FP_TYPE* x, int nx) {
  return nx % 2 == 1 ? cig_qselect(x, nx, nx / 2) :
    (cig_qselect(x, nx, nx / 2 - 1) + cig_qselect(x, nx, nx / 2)) * 0.5;
}

FP_TYPE* cig_xcorr(FP_TYPE* x, FP_TYPE* y, int nx, int maxlag) {
  FP_TYPE* R = calloc(maxlag, sizeof(FP_TYPE));
  for(int m = 0; m < maxlag; m ++) {
    for(int i = 0; i < nx - m; i ++)
      R[m] += x[m + i] * y[i];
  }
  return R;
}

// orient: 1 (maximum) or -1 (minimum)
int cig_find_peak(FP_TYPE* x, int lidx, int uidx, int orient) {
  FP_TYPE max = x[lidx] * orient;
  FP_TYPE max_idx = lidx;
  for(int i = lidx; i <= uidx; i ++)
    if(x[i] * orient > x[i - 1] * orient && x[i] * orient > x[i + 1] * orient)
      if(x[i] * orient > max) {
        max = x[i] * orient;
        max_idx = i;
      }
  return max_idx;
}

FP_TYPE* cig_gensins(FP_TYPE* freq, FP_TYPE* ampl, FP_TYPE* phse,
  int nsin, int fs, int n) {
  FP_TYPE* x = calloc(n, sizeof(FP_TYPE));
  FP_TYPE* s = calloc(n, sizeof(FP_TYPE));
  for(int i = 0; i < nsin; i ++) {
    FP_TYPE i_f = freq[i];
    FP_TYPE i_a = ampl[i];
    FP_TYPE i_p = phse[i];

    FP_TYPE tpffs = 2.0 * M_PI / fs * i_f;
    FP_TYPE c = 2.0 * cos_3(tpffs);
    s[0] = cos_3(tpffs * (- n / 2) + i_p) * i_a;
    s[1] = cos_3(tpffs * (- n / 2 + 1) + i_p) * i_a;
    x[0] += s[0];
    x[1] += s[1];
    for(int t = 2; t < n; t ++) {
      s[t] = c * s[t - 1] - s[t - 2];
      x[t] += s[t];
    }
  }
  free(s);
  return x;
}

static FP_TYPE* get_window(char* name, int nw, int optlv) {
  if(optlv >= 3) {
    if(! strcmp(name, "hanning"))
      return hanning(nw);
    else if(! strcmp(name, "hamming"))
      return hamming(nw);
    else if(! strcmp(name, "blackman"))
      return blackman(nw);
    else if(! strcmp(name, "blackman_harris"))
      return blackman_harris(nw);
    return boxcar(nw);
  } else {
    if(! strcmp(name, "hanning"))
      return hanning_2(nw);
    else if(! strcmp(name, "hamming"))
      return hamming_2(nw);
    else if(! strcmp(name, "blackman"))
      return blackman_2(nw);
    else if(! strcmp(name, "blackman_harris"))
      return blackman_harris_2(nw);
    return boxcar(nw);
  }
}

FP_TYPE* cig_dct(FP_TYPE* x, int nx) {
  FP_TYPE* ret = calloc(nx, sizeof(FP_TYPE));
  for(int k = 0; k < nx; k ++)
    for(int n = 0; n < nx; n ++)
      ret[k] += x[n] * cos_2(M_PI / nx * (n + 0.5) * k);
  ret[0] /= sqrt(nx);
  for(int k = 1; k < nx; k ++)
    ret[k] /= sqrt(nx * 0.5);
  return ret;
}

void cig_fft(FP_TYPE* xr, FP_TYPE* xi, FP_TYPE* yr, FP_TYPE* yi,
    int n, FP_TYPE* buffer, FP_TYPE mode) {
  for(int i = 0; i < n; i ++) {
    buffer[i * 2] = xr == NULL ? 0 : xr[i];
    buffer[i * 2 + 1] = xi == NULL ? 0 : xi[i];
  }
  cdft(2 * n, mode, buffer);
  if(mode < 0)
    for(int i = 0; i < n; i ++) {
      if(yr != NULL) yr[i] = buffer[i * 2];
      if(yi != NULL) yi[i] = buffer[i * 2 + 1];
    }
  else
    for(int i = 0; i < n; i ++) {
      if(yr != NULL) yr[i] = buffer[i * 2] / n;
      if(yi != NULL) yi[i] = buffer[i * 2 + 1] / n;
    }
}

void cig_idft(FP_TYPE* xr, FP_TYPE* xi, FP_TYPE* yr, FP_TYPE* yi, int n) {
  FP_TYPE tpon = 2.0 * M_PI / n;
  for(int k = 0; k < n; k ++) {
    FP_TYPE sumr = 0;
    FP_TYPE sumi = 0;
    FP_TYPE tponk = tpon * k;
    if(xr != NULL && xi != NULL)
      for(int i = 0; i < n; i ++) {
        FP_TYPE re = cos_3(tponk * i);
        FP_TYPE im = sin_3(tponk * i);
        sumr += xr[i] * re - xi[i] * im;
        sumi += xr[i] * im + xi[i] * re;
      }
    else if(xr != NULL && xi == NULL)
      for(int i = 0; i < n; i ++) {
        FP_TYPE re = cos_3(tponk * i);
        FP_TYPE im = sin_3(tponk * i);
        sumr += xr[i] * re;
        sumi += xr[i] * im;
      }
    else if(xr == NULL && xi != NULL)
      for(int i = 0; i < n; i ++) {
        FP_TYPE re = cos_3(tponk * i);
        FP_TYPE im = sin_3(tponk * i);
        sumr -= xi[i] * im;
        sumi += xi[i] * re;
      }
    if(yr != NULL) yr[k] = sumr / n;
    if(yi != NULL) yi[k] = sumi / n;
  }
}

FP_TYPE* cig_levinson(FP_TYPE* R, int n) {
  FP_TYPE* a = calloc(n, sizeof(FP_TYPE));
  FP_TYPE* a_new = calloc(n, sizeof(FP_TYPE));
  FP_TYPE k = -R[1] / (R[0] * 1.000001);
  a[0] = 1.0;
  a[1] = k;
  FP_TYPE gain = R[0] * 1.000001 * (1 - k * k);
  for(int i = 2; i < n; i ++) {
    FP_TYPE s = R[i];
    for(int j = 1; j < i; j ++)
      s += R[j] * a[i - j];
    k = -s / gain;
    for(int j = 1; j < i; j ++)
      a_new[j] = a[j] + k * a[i - j];
    for(int j = 1; j < i; j ++)
      a[j] = a_new[j];
    a[i] = k;
    gain *= 1.0 - k * k;
  }
  free(a_new);
  return a;
}

FP_TYPE* cig_winfir(int order, FP_TYPE cutoff, FP_TYPE cutoff2, char* type, char* window) {
  FP_TYPE cutk  = cutoff  * order;
  FP_TYPE cutk2 = cutoff2 * order;
  FP_TYPE* freqrsp = calloc(order, sizeof(FP_TYPE));
  FP_TYPE* timersp = calloc(order, sizeof(FP_TYPE));

  FP_TYPE* w = get_window(window, order, 3);
  
  if(! strcmp(type, "lowpass")) {
    for(int i = 0; i <= floor(cutk); i ++)
      freqrsp[i] = 1;
    freqrsp[(int)ceil(cutk)] = fmod(cutk, 1.0);
  } else
  if(! strcmp(type, "highpass")) {
    for(int i = order / 2; i > floor(cutk); i --)
      freqrsp[i] = 1;
    freqrsp[(int)floor(cutk)] = 1.0 - fmod(cutk, 1.0);
  } else
  if(! strcmp(type, "bandpass")) {
    for(int i = floor(cutk2); i > floor(cutk); i --)
      freqrsp[i] = 1;
    freqrsp[(int)floor(cutk)] = 1.0 - fmod(cutk, 1.0);
    freqrsp[(int)ceil(cutk2)] = fmod(cutk2, 1.0);
  } else
  if(! strcmp(type, "bandstop")) {
    for(int i = 0; i < order; i ++) freqrsp[i] = 1;
    for(int i = floor(cutk2); i > floor(cutk); i --)
      freqrsp[i] = 0;
    freqrsp[(int)floor(cutk)] = fmod(cutk, 1.0);
    freqrsp[(int)ceil(cutk2)] = 1.0 - fmod(cutk2, 1.0);
  }
  
  complete_symm(freqrsp, order);
  idft(freqrsp, NULL, timersp, NULL, order);
  FP_TYPE* h = fftshift(timersp, order);
  for(int i = 0; i < order; i ++)
    h[i] *= w[i];
  
  free(w);
  free(freqrsp);
  free(timersp);
  return h;
}

FP_TYPE* cig_convolution(FP_TYPE* x, FP_TYPE* h, int nx, int nh) {
  FP_TYPE* xpad = calloc(nx + nh * 2 - 1, sizeof(FP_TYPE));
  FP_TYPE* y = calloc(nx + nh - 1, sizeof(FP_TYPE));
  for(int i = 0; i < nx; i ++)
    xpad[i + nh - 1] = x[i];
  for(int i = 0; i < nx + nh - 1; i ++)
    for(int k = 0; k < nh; k ++)
      y[i] += h[k] * xpad[i + nh - 1 - k];
  free(xpad);
  return y;
}

// unrolled filters

static FP_TYPE* cig_filter_order6(FP_TYPE* b, FP_TYPE* a, FP_TYPE* x, int nx) {
  FP_TYPE* y = calloc(nx + 5, sizeof(FP_TYPE));
  for(int i = 5; i < nx; i ++) {
      y[i] -= a[1] * y[i - 1] + a[2] * y[i - 2] + a[3] * y[i - 3] + a[4] * y[i - 4] + a[5] * y[i - 5];
      y[i] += b[0] * x[i - 0] + b[1] * x[i - 1] + b[2] * x[i - 2] + b[3] * x[i - 3] + b[4] * x[i - 4] + b[5] * x[i - 5];
  }
  return y;
}

static FP_TYPE* cig_filter_order5(FP_TYPE* b, FP_TYPE* a, FP_TYPE* x, int nx) {
  FP_TYPE* y = calloc(nx + 4, sizeof(FP_TYPE));
  for(int i = 4; i < nx; i ++) {
      y[i] -= a[1] * y[i - 1] + a[2] * y[i - 2] + a[3] * y[i - 3] + a[4] * y[i - 4];
      y[i] += b[0] * x[i - 0] + b[1] * x[i - 1] + b[2] * x[i - 2] + b[3] * x[i - 3] + b[4] * x[i - 4];
  }
  return y;
}

static FP_TYPE* cig_filter_order4(FP_TYPE* b, FP_TYPE* a, FP_TYPE* x, int nx) {
  FP_TYPE* y = calloc(nx + 3, sizeof(FP_TYPE));
  for(int i = 3; i < nx; i ++) {
      y[i] -= a[1] * y[i - 1] + a[2] * y[i - 2] + a[3] * y[i - 3];
      y[i] += b[0] * x[i - 0] + b[1] * x[i - 1] + b[2] * x[i - 2] + b[3] * x[i - 3];
  }
  return y;
}

static FP_TYPE* cig_filter_order3(FP_TYPE* b, FP_TYPE* a, FP_TYPE* x, int nx) {
  FP_TYPE* y = calloc(nx + 2, sizeof(FP_TYPE));
  for(int i = 2; i < nx; i ++) {
      y[i] -= a[1] * y[i - 1] + a[2] * y[i - 2];
      y[i] += b[0] * x[i - 0] + b[1] * x[i - 1] + b[2] * x[i - 2];
  }
  return y;
}

static FP_TYPE* cig_filter_order2(FP_TYPE* b, FP_TYPE* a, FP_TYPE* x, int nx) {
  FP_TYPE* y = calloc(nx + 1, sizeof(FP_TYPE));
  for(int i = 1; i < nx; i ++) {
      y[i] -= a[1] * y[i - 1];
      y[i] += b[0] * x[i - 0] + b[1] * x[i - 1];
  }
  return y;
}

FP_TYPE* cig_filter(FP_TYPE* b, int nb, FP_TYPE* a, int na, FP_TYPE* x, int nx) {
  if(na == nb && na < 7) {
    if(na == 6)
      return cig_filter_order6(b, a, x, nx);
    if(na == 5)
      return cig_filter_order5(b, a, x, nx);
    if(na == 4)
      return cig_filter_order4(b, a, x, nx);
    if(na == 3)
      return cig_filter_order3(b, a, x, nx);
    if(na == 2)
      return cig_filter_order2(b, a, x, nx);
  }

  int nh = max(na, nb);
  FP_TYPE* y = calloc(nx + nh - 1, sizeof(FP_TYPE));
  for(int i = 0; i < nh; i ++) {
    for(int k = 1; k < na; k ++)
      if(i - k >= 0)
        y[i] -= a[k] * y[i - k];
    for(int k = 0; k < nb; k ++)
      if(i - k >= 0)
        y[i] += b[k] * x[i - k];
  }
  for(int i = nh; i < nx; i ++) {
    for(int k = 1; k < na; k ++)
      y[i] -= a[k] * y[i - k];
    for(int k = 0; k < nb; k ++)
      y[i] += b[k] * x[i - k];
  }
  return y;
}

static FP_TYPE interp_kernel(FP_TYPE x, FP_TYPE a) {
  if(x == 0) return 1;
  if(x > a || x < -a) return 0;
  return a * sin_2(M_PI * x) * sin_2(M_PI * x / a) / M_PI / M_PI / x / x;
}

// xi should be ascending while the choice of x can be arbitrary
FP_TYPE* cig_interp(FP_TYPE* xi, FP_TYPE* yi, int ni, FP_TYPE* x, int nx) {
  FP_TYPE* y = malloc(nx * sizeof(FP_TYPE));
  int srcidx = 0;
  for(int i = 0; i < nx; i ++) {
    FP_TYPE dstx = x[i];
    while(srcidx + 1 < ni && xi[srcidx + 1] < dstx) srcidx ++;
    int i0 = srcidx == 0 ? 0 : srcidx;
    int i1 = srcidx == ni - 1 ? srcidx : srcidx + 1;
    if(i0 != i1 && dstx > xi[0])
      y[i] = (yi[i1] - yi[i0]) * (dstx - xi[i0]) / (xi[i1] - xi[i0]) + yi[i0];
    else
      y[i] = yi[i0];
  }
  return y;
}

// interpolation on a uniformly sampled signal, assuming x is ascending
FP_TYPE* cig_interpu(FP_TYPE xi0, FP_TYPE xi1, FP_TYPE* yi, int ni, FP_TYPE* x, int nx) {
  FP_TYPE* y = malloc(nx * sizeof(FP_TYPE));
  int begin = 0, end = nx - 1;
  while(x[begin] < xi0 && begin < nx) {
    y[begin] = yi[0];
    begin ++;
  }
  for(end = nx - 1; end > begin; end --) {
    FP_TYPE srcidx = (x[end] - xi0) / (xi1 - xi0) * ni;
    if(srcidx + 1.0 < ni) break;
    y[end] = yi[ni - 1];
  }
  for(int i = begin; i <= end; i ++) {
    FP_TYPE srcidx = (x[i] - xi0) / (xi1 - xi0) * ni;
    int base = srcidx;
    FP_TYPE r = srcidx - base;
    y[i] = yi[base] + (yi[base + 1] - yi[base]) * r;
  }
  return y;
}

FP_TYPE* cig_sincinterpu(FP_TYPE xi0, FP_TYPE xi1, FP_TYPE* yi, int ni, FP_TYPE* x, int nx) {
  FP_TYPE* y = malloc(nx * sizeof(FP_TYPE));
  for(int i = 0; i < nx; i ++) {
    FP_TYPE srcidx = (x[i] - xi0) / (xi1 - xi0) * ni;
    int ix_base = floor(srcidx);
    FP_TYPE a = srcidx - ix_base;
    for(int j = -3; j <= 4; j ++)
      if(ix_base + j > 0 && ix_base + j < ni)
        y[i] += interp_kernel(a - j, 4) * yi[ix_base + j];
  }
  return y;
}

typedef FP_TYPE (*fp_n_to_one)(FP_TYPE* x);

FP_TYPE* cig_medfilt(FP_TYPE* x, int nx, int order) {
  fp_n_to_one ffilt = NULL;
  if(order == 3)
    ffilt = opt_med3;
  else if(order == 5)
    ffilt = opt_med5;
  else if(order == 7)
    ffilt = opt_med7;
  else if(order == 9)
    ffilt = opt_med9;
  else if(order == 25)
    ffilt = opt_med25;

  FP_TYPE* y = malloc(nx * sizeof(FP_TYPE));
  FP_TYPE* tmp = malloc(order * sizeof(FP_TYPE));
  int halford = order / 2;
  for(int i = 0; i < halford; i ++)
    y[i] = medianfp(x, i + 1);
  for(int i = 0; i < order; i ++)
    tmp[i] = x[i];
  for(int i = halford; i < nx - halford; i ++) {
    for(int k = 0; k < order; k ++)
      tmp[k] = x[i + k - halford];
    if(ffilt != NULL)
      y[i] = ffilt(tmp);
    else
      y[i] = medianfp(tmp, order);
  }
  for(int i = nx - halford; i < nx; i ++)
    y[i] = medianfp(x + i, nx - i);
  free(tmp);
  return y;
}

FP_TYPE* cig_rresample(FP_TYPE* x, int nx, FP_TYPE ratio, int* ny) {
  *ny = ratio == 1.0 ? nx : round(nx * ratio);
  FP_TYPE* xlp = x;
  FP_TYPE* y = calloc(*ny, sizeof(FP_TYPE));
  if(ratio == 1.0) {
    memcpy(y, x, sizeof(FP_TYPE) * nx);
    return y;
  }

  if(ratio < 1.0) { // low pass the signal first
    FP_TYPE* h = fir1(32, ratio, "lowpass", "hamming");
    xlp = conv(x, h, nx, 32) + 16;
    free(h);
  }

  for(int i = 0; i < *ny; i ++) {
    FP_TYPE ix = (FP_TYPE)i / ratio;
    int ix_base = floor(ix);
    FP_TYPE a = ix - ix_base;
    for(int j = -3; j <= 4; j ++)
      if(ix_base + j > 0 && ix_base + j < nx)
        y[i] += interp_kernel(a - j, 4) * xlp[ix_base + j];
  }

  if(ratio < 1.0) {
    xlp -= 16;
    free(xlp);
  }

  return y;
}

static int sign(FP_TYPE x) {
  return x > 0.0 ? 1 : (x == 0.0 ? 0 : -1);
}

static FP_TYPE fzero_(fpe_one_to_one func, FP_TYPE xmin, FP_TYPE xmax, int depth, void* env) {
  FP_TYPE xmean = (xmin + xmax) * 0.5;
  FP_TYPE ymean = func(xmean, env);
  if(fabs(ymean) < 1e-5 || depth > 100) return xmean;
  FP_TYPE ymax = func(xmax, env);
  if(sign(ymax) != sign(ymean)) return fzero_(func, xmean, xmax, depth + 1, env);
  return fzero_(func, xmin, xmean, depth + 1, env);
}

FP_TYPE cig_fzero(fpe_one_to_one func, FP_TYPE xmin, FP_TYPE xmax, void* env) {
  if(sign(func(xmax, env)) == sign(func(xmin, env))) return 0;
  return fzero_(func, xmin, xmax, 0, env);
}

// Pearson Correlation Coefficient
static FP_TYPE corr_kernel_acf(
  FP_TYPE* x, FP_TYPE* xx, double* x1, double* x2, int w, int d) {
  FP_TYPE a = 0;
  FP_TYPE b = 0;
  FP_TYPE sumx = x1[w] - x1[0];
  FP_TYPE sumxd = x1[w + d] - x1[d];
  FP_TYPE sumx2 = x2[w] - x2[0];
  FP_TYPE sumx2d = x2[w + d] - x2[d];
  for(int i = 0; i < w; i ++)
    a += x[i] * x[i + d];
  a *= w;
  a -= sumx * sumxd;
  b = (w * sumx2 - sumx * sumx) * (w * sumx2d - sumxd * sumxd);
  FP_TYPE r = a / sqrt(b);
  r = max(r, 0); r = min(r, 1);
  return r;
}

static FP_TYPE corr_kernel_amdf(
  FP_TYPE* x, FP_TYPE* xx, double* x1, double* x2, int w, int d) {
  FP_TYPE a = 0;
  for(int i = 0; i < w; i ++)
    a += fabs(x[i] - x[i + d]);
  return a;
}

void cig_correlogram(FP_TYPE* x, int nx, int* center, int* nwin, int nfrm,
  int max_period, int method, FP_TYPE** R) {
  double* x1 = calloc(nx + 1, sizeof(double)); // integrate x[i]
  double* x2 = calloc(nx + 1, sizeof(double)); // integrate x^2[i]
  FP_TYPE* xx = calloc(nx, sizeof(FP_TYPE));   // x^2[i]
  for(int i = 0; i < nx; i ++)
    x1[i + 1] = x1[i] + (double)x[i];
  for(int i = 0; i < nx; i ++)
    x2[i + 1] = x2[i] + (double)x[i] * (double)x[i];
  for(int i = 0; i < nx; i ++)
    xx[i] = x[i] * x[i];
  for(int i = 0; i < nfrm; i ++) {
    int t = center[i];
    int w = nwin[i];
    t = max(0, t);
    t = min(nx - 1, t);
    for(int d = 0; d < max_period; d ++)
      if(t + w + d < nx) {
        if(method == CIG_CORR_ACF)
          R[i][d] = corr_kernel_acf(x + t, xx + t, x1 + t, x2 + t, w, d);
        else if(method == CIG_CORR_AMDF)
          R[i][d] = corr_kernel_amdf(x + t, xx + t, x1 + t, x2 + t, w, d);
      } else
        R[i][d] = 0;
  }
  free(x1);
  free(x2);
  free(xx);
}

FP_TYPE** cig_invcrgm(FP_TYPE** R, int nfrm, int max_period, int fs, FP_TYPE* faxis, int nf) {
  FP_TYPE** Ri = malloc2d(nfrm, nf, sizeof(FP_TYPE*));
  FP_TYPE* Raxis = linspace(0, max_period - 1, max_period);
  FP_TYPE* invmap = calloc(nf, sizeof(FP_TYPE));
  for(int i = 0; i < nf; i ++)
    invmap[i] = fs / faxis[nf - i - 1];
  for(int i = 0; i < nfrm; i ++) {
    FP_TYPE* rflip = interp1(Raxis, R[i], max_period, invmap, nf);
    for(int j = 0; j < nf; j ++)
      Ri[i][j] = rflip[nf - j - 1];
    free(rflip);
  }
  free(invmap);
  free(Raxis);
  return Ri;
}

void cig_stft_forward(FP_TYPE* x, int nx, int* center, int* nwin, int nfrm,
  int nfft, char* window, int subt_mean, int optlv,
  FP_TYPE* norm_factor, FP_TYPE* weight_factor, FP_TYPE** Xmagn, FP_TYPE** Xphse) {

# ifndef _OPENMP
  FP_TYPE* buff = calloc(nfft * 5, sizeof(FP_TYPE));
  FP_TYPE* fftbuff = buff;
  FP_TYPE* xbuff = buff + nfft * 2;
  FP_TYPE* ybuffr = buff + nfft * 3;
  FP_TYPE* ybuffi = buff + nfft * 4;
# endif

  FP_TYPE* w = NULL;
  if(norm_factor != NULL || weight_factor != NULL)
    w = get_window(window, nwin[0], optlv);

  if(weight_factor != NULL) {
    *weight_factor = 0;
    for(int i = 0; i < nwin[0]; i ++)
      *weight_factor += w[i];
  }
  if(norm_factor != NULL) {
    *norm_factor = 0;
    for(int i = 0; i < nwin[0]; i += center[1] - center[0])
      *norm_factor += w[i];
  }

# ifdef _OPENMP
# pragma omp parallel for
# endif
  for(int t = 0; t < nfrm; t ++) {
#   ifdef _OPENMP
    FP_TYPE* buff = calloc(nfft * 5, sizeof(FP_TYPE));
    FP_TYPE* fftbuff = buff;
    FP_TYPE* xbuff = buff + nfft * 2;
    FP_TYPE* ybuffr = buff + nfft * 3;
    FP_TYPE* ybuffi = buff + nfft * 4;
#   endif
    int tn = center[t];
    
    FP_TYPE* xfrm;
    FP_TYPE* spec_magn, *spec_phse;
    xfrm = spec_magn = spec_phse = NULL;
    
    FP_TYPE* wlocal = w;
    if(norm_factor == NULL && weight_factor == NULL) {
      wlocal = get_window(window, nwin[t], optlv);
    }

    xfrm = fetch_frame(x, nx, tn, nwin[t]);
    if(subt_mean) {
      FP_TYPE mean_xfrm = sumfp(xfrm, nwin[t]) / nwin[t];
      for(int i = 0; i < nwin[t]; i ++)
        xfrm[i] = (xfrm[i] - mean_xfrm) * wlocal[i];
    } else {
      for(int i = 0; i < nwin[t]; i ++)
        xfrm[i] *= wlocal[i];
    }

    memset(xbuff, 0, nfft * sizeof(FP_TYPE));
    for(int i = 0; i < nwin[t] / 2; i ++) {
      xbuff[i] = xfrm[i + nwin[t] / 2];
      xbuff[nfft - nwin[t] / 2 + i] = xfrm[i];
    }
    fft(xbuff, NULL, ybuffr, ybuffi, nfft, fftbuff);
    
    if(Xmagn != NULL) {
      spec_magn = abscplx(ybuffr, ybuffi, nfft / 2 + 1);
      for(int i = 0; i < nfft / 2 + 1; i ++)
        Xmagn[t][i] = spec_magn[i];
      free(spec_magn);
    }

    if(Xphse != NULL) {
      spec_phse = argcplx(ybuffr, ybuffi, nfft / 2 + 1);
      for(int i = 0; i < nfft / 2 + 1; i ++)
        Xphse[t][i] = spec_phse[i];
      free(spec_phse);
    }

    free(xfrm);
    if(norm_factor == NULL && weight_factor == NULL)
      free(wlocal);

#   ifdef _OPENMP
    free(buff);
#   endif
  }
  
  if(norm_factor != NULL || weight_factor != NULL)
    free(w);
# ifndef _OPENMP
  free(buff);
# endif
}

FP_TYPE* cig_stft_backward(FP_TYPE** Xmagn, FP_TYPE** Xphse, int nhop, int nfrm,
  int offset, int hop_factor, int zp_factor, int nfade, FP_TYPE norm_factor, int* ny) {
  int nwin = nhop * hop_factor;
  int nfft = nwin * zp_factor;
  FP_TYPE* wfade = hanning(nfade * 2);
  
  *ny = nhop * nfrm + offset;
  FP_TYPE* y = calloc(*ny, sizeof(FP_TYPE));
  
# ifndef _OPENMP
  FP_TYPE* buff = calloc(nfft * 5, sizeof(FP_TYPE));
  FP_TYPE* fftbuff = buff;
  FP_TYPE* ybuff = buff + nfft * 2;
  FP_TYPE* xbuffr = buff + nfft * 3;
  FP_TYPE* xbuffi = buff + nfft * 4;
# endif
  
# ifdef _OPENMP
# pragma omp parallel for
# endif
  for(int t = 0; t < nfrm; t ++) {
#   ifdef _OPENMP
    FP_TYPE* buff = calloc(nfft * 5, sizeof(FP_TYPE));
    FP_TYPE* fftbuff = buff;
    FP_TYPE* ybuff = buff + nfft * 2;
    FP_TYPE* xbuffr = buff + nfft * 3;
    FP_TYPE* xbuffi = buff + nfft * 4;
#   endif

    int tn = t * nhop + offset;
    
    for(int i = 0; i < nfft / 2 + 1; i ++) {
      xbuffr[i] = Xmagn[t][i] * cos_2(Xphse[t][i]);
      xbuffi[i] = Xmagn[t][i] * sin_2(Xphse[t][i]);
    }
    complete_symm (xbuffr, nfft);
    complete_asymm(xbuffi, nfft);
    ifft(xbuffr, xbuffi, ybuff, NULL, nfft, fftbuff);
    
    FP_TYPE* yfrm = fftshift(ybuff, nfft);
    for(int i = 0; i < nfade; i ++) {
      yfrm[i] *= wfade[i];
      yfrm[nfft - i - 1] *= wfade[i];
    }
    
#   ifdef _OPENMP
#   pragma omp critical
#   endif
    for(int i = 0; i < nfft; i ++) {
      int idx = tn + i - nfft / 2;
      if(idx >= 0 && idx < *ny) {
        y[idx] += yfrm[i] / norm_factor;
      }
    }
    
    free(yfrm);
#   ifdef _OPENMP
    free(buff);
#   endif
  }
  
  free(wfade);
# ifndef _OPENMP
  free(buff);
# endif
  return y;
}

FP_TYPE cig_qifft(FP_TYPE* magn, int k, FP_TYPE* dst_freq) {
  FP_TYPE a, b, c, a1, a2, x;
  a = magn[k - 1];
  b = magn[k + 0];
  c = magn[k + 1];
  a1 = (a + c) / 2.0 - b;
  a2 = c - b - a1;
  x = - a2 / a1 * 0.5;
  
  x = (fabs(x) < 1.0) ? x : 0; // in case we get some x outside of [k-1, k+1]
  
  *dst_freq = (FP_TYPE)k + x;
  FP_TYPE ret = a1 * x * x + a2 * x + b;
  return ret > b + 0.2 ? b + 0.2 : ret; // in case it pops
}

static inline FP_TYPE barkmask(FP_TYPE z) {
  if(z < -1.3 || z > 2.5) return 0;
  if(z < -0.5) return fastpow(10.0, 2.5 * (z + 0.5));
  if(z <= 0.5) return 1.0;
  return fastpow(10.0, (0.5 - z));
}

filterbank* cig_create_empty_filterbank(int nf, FP_TYPE fnyq, int nchannel) {
  filterbank* ret = malloc(sizeof(filterbank));
  ret -> nchannel = nchannel;
  ret -> nf = nf;
  ret -> fnyq = fnyq;
  ret -> fresp = calloc(nchannel, sizeof(FP_TYPE*));
  ret -> lower_idx = calloc(nchannel, sizeof(int));
  ret -> upper_idx = calloc(nchannel, sizeof(int));
  for(int i = 0; i < nchannel; i ++) {
    ret -> fresp[i] = calloc(nf, sizeof(FP_TYPE));
    ret -> upper_idx[i] = nf;
  }
  return ret;
}

filterbank* cig_create_plp_filterbank(int nf, FP_TYPE fnyq, int nchannel) {
  filterbank* ret = cig_create_empty_filterbank(nf, fnyq, nchannel);
  for(int i = 0; i < nchannel; i ++) {
    FP_TYPE weight = eqloud(bark2freq(i));
    for(int k = 0; k < nf; k ++)
      ret -> fresp[i][k] = barkmask(i - freq2bark((FP_TYPE)k / nf * fnyq)) * weight;
  }
  return ret;
}

filterbank* cig_create_melfreq_filterbank(int nf, FP_TYPE fnyq, int nchannel,
  FP_TYPE min_freq, FP_TYPE max_freq, FP_TYPE scale) {
  filterbank* ret = cig_create_empty_filterbank(nf, fnyq, nchannel);
  
  FP_TYPE* freqs = melspace(min_freq, max_freq, nchannel);
  for(int j = 0; j < nchannel; j ++) {
    FP_TYPE f_0 = j == 0 ? 0 : freqs[j - 1];
    FP_TYPE f_1 = freqs[j];
    FP_TYPE f_2 = freqs[j + 1];
    if(f_0 > f_1 - 400)
      f_0 = max(0, f_1 - 400);
    if(f_2 < f_1 + 400)
      f_2 = f_1 + 400;
    int lower_idx = j == 0 ? 0 : floor(f_0 * nf / fnyq * scale);
    int upper_idx = ceil(f_1 * nf / fnyq * scale);
    int centr_idx = round(f_2 * nf / fnyq * scale);
    
    ret -> lower_idx[j] = lower_idx;
    ret -> upper_idx[j] = upper_idx;
    for(int k = lower_idx; k < centr_idx; k ++)
      ret -> fresp[j][k] = (FP_TYPE)(k - lower_idx + 1) / (centr_idx - lower_idx);
    for(int k = centr_idx; k < upper_idx; k ++)
      ret -> fresp[j][k] = (1.0 - (FP_TYPE)(k - centr_idx) / (upper_idx - centr_idx));
  }
  free(freqs);
  return ret;
}

void cig_delete_filterbank(filterbank* dst) {
  if(dst == NULL) return;
  for(int i = 0; i < dst -> nchannel; i ++)
    free(dst -> fresp[i]);
  free(dst -> lower_idx);
  free(dst -> upper_idx);
  free(dst -> fresp);
  free(dst);
}

FP_TYPE** cig_filterbank_spectrogram(filterbank* fbank, FP_TYPE** S, int nfrm,
  int nfft, int fs, int crtenergy) {
  FP_TYPE** X = malloc2d(nfrm, fbank -> nchannel, sizeof(FP_TYPE));
# ifdef _OPENMP
# pragma omp parallel for
# endif
  for(int i = 0; i < nfrm; i ++) {
    for(int j = 0; j < fbank -> nchannel; j ++) {
      X[i][j] = 0;
      for(int k = fbank -> lower_idx[j]; k < fbank -> upper_idx[j]; k ++)
        X[i][j] += fbank -> fresp[j][k] * S[i][k];
      X[i][j] = crtenergy ? pow(X[i][j], 0.33) : log_2(X[i][j] + M_EPS);
    }
  }
  return X;
}

FP_TYPE* cig_filterbank_spectrum(filterbank* fbank, FP_TYPE* S, int nfft, int fs,
  int crtenergy) {
  FP_TYPE* X = calloc(fbank -> nchannel, sizeof(FP_TYPE));
  for(int j = 0; j < fbank -> nchannel; j ++) {
    for(int k = fbank -> lower_idx[j]; k < fbank -> upper_idx[j]; k ++)
      X[j] += fbank -> fresp[j][k] * S[k];
    X[j] = crtenergy ? pow(X[j], 0.33) : log_2(X[j] + M_EPS);
  }
  return X;
}

// Morise, Masanori. "Cheaptrick, a spectral envelope estimator for high-quality
//   speech synthesis." Speech Communication 67 (2015): 1-7.
FP_TYPE* cig_spec2env(FP_TYPE* S, int nfft, int fs, FP_TYPE f0, FP_TYPE* Cout) {
  FP_TYPE* buff = malloc(nfft * 4 * sizeof(FP_TYPE));
  FP_TYPE* V = buff;
  FP_TYPE* C = buff + nfft;
  FP_TYPE* fftbuff = buff + nfft * 2;
  FP_TYPE smoothord = (FP_TYPE)nfft / (fs / f0 * 3.0);

  for(int i = 0; i < nfft / 2 + 1; i ++) V[i] = S[i];
  int kf0 = ceil(f0 / fs * nfft);
  for(int i = 0; i < kf0 / 2; i ++) V[i] = V[kf0 - i];
  FP_TYPE* smoothed = moving_avg(V, nfft / 2 + 1, smoothord);

  for(int i = 0; i < nfft / 2 + 1; i ++)
    V[i] = log_2(smoothed[i] + M_EPS);
  complete_symm(V, nfft);
  ifft(V, NULL, C, NULL, nfft, fftbuff);
  for(int i = 1; i < nfft / 2 + 1; i ++)
    C[i] *= sin_2(i * f0 / fs * M_PI) / (i * f0 / fs * M_PI) *
      (1.18 - 2.0 * 0.09 * cos_2(2.0 * M_PI * i * f0 / fs));
  complete_symm(C, nfft);
  fft(C, NULL, V, NULL, nfft, fftbuff);
  free(smoothed);
  if(Cout != NULL)
    for(int i = 0; i < nfft / 2 + 1; i ++)
      Cout[i] = C[i];
  return realloc(V, nfft * sizeof(FP_TYPE));
}

// Huber, Stefan, and Axel Roebel. "On the use of voice descriptors for glottal
//   source shape parameter estimation." Computer Speech & Language 28.5 (2014):
//   1170-1194.
// Fant, Gunnar. "The LF-model revisited. Transformations and frequency domain
//    analysis." Speech Trans. Lab. Q. Rep., Royal Inst. of Tech. Stockholm 2.3
//    (1995): 40.
lfmodel cig_lfmodel_from_rd(FP_TYPE rd, FP_TYPE T0, FP_TYPE Ee) {
  lfmodel ret;
  FP_TYPE Rap = rd < 0.21 ? 1e-6 : (rd < 2.7 ? (-1.0 + 4.8 * rd) / 100.0 : 0.323 / rd);
  FP_TYPE OQupp = 1.0 - 1.0 / (2.17 * rd);
  FP_TYPE Rkp, Rgp;
  if(rd < 2.7) {
    Rkp = (22.4 + 11.8 * rd) / 100.0;
    Rgp = 0.25 * Rkp / ((0.11 * rd) / (0.5 + 1.2 * Rkp) - Rap);
  } else {
    Rgp = 9.3552e-3 + 596e-2 / (7.96 - 2.0 * OQupp);
    Rkp = 2.0 * Rgp * OQupp - 1.0428;
  }
  ret.tp = 1.0 / (2.0 * Rgp);
  ret.te = ret.tp * (Rkp + 1.0);
  ret.ta = Rap;
  ret.T0 = T0;
  ret.Ee = Ee;
  return ret;
};

typedef struct {
  FP_TYPE T0;
  FP_TYPE Te;
  FP_TYPE Tp;
  FP_TYPE Ta;
  FP_TYPE wg;
  FP_TYPE sin_wgTe;
  FP_TYPE cos_wgTe;
  FP_TYPE e;
  FP_TYPE A;
  FP_TYPE a;
  FP_TYPE E0;
} lfparam;

static FP_TYPE efunc(FP_TYPE x, void* env) {
  lfparam* tmpparam = (lfparam*)env;
  return 1.0 - exp_3((tmpparam -> Te - tmpparam -> T0) * x) - tmpparam -> Ta * x;
}

static FP_TYPE afunc(FP_TYPE x, void* env) {
  lfparam* tmpparam = (lfparam*)env;
  return (x * x + tmpparam -> wg * tmpparam -> wg) * tmpparam -> sin_wgTe * tmpparam -> A +
    tmpparam -> wg * exp_3(- x * tmpparam -> Te) + x * tmpparam -> sin_wgTe -
    tmpparam -> wg * tmpparam -> cos_wgTe;
}

// Fant, Gunnar, Johan Liljencrants, and Qi-guang Lin. "A four-parameter model
//    of glottal flow." STL-QPSR 4.1985 (1985): 1-13.
static lfparam lfparam_from_lfmodel(lfmodel model) {
  lfparam ret = {
    .T0 = model.T0,
    .Te = model.T0 * model.te,
    .Tp = model.T0 * model.tp,
    .Ta = model.T0 * model.ta
  };
  ret.wg = M_PI / ret.Tp;
  ret.sin_wgTe = sin_3(ret.wg * ret.Te);
  ret.cos_wgTe = cos_3(ret.wg * ret.Te);
  FP_TYPE e = fzero(efunc, 1.0, 2.0 / (ret.Ta + 1e-9), & ret);
  FP_TYPE e_Te_T0 = exp_3(e * (ret.Te - ret.T0));
  ret.A = (1.0 - e_Te_T0) / (e * e * ret.Ta) + (ret.Te - ret.T0) * e_Te_T0 / (e * ret.Ta);
  ret.e = e;
  ret.a = fzero(afunc, 0.0, 1e9, & ret);
  ret.E0 = -model.Ee / exp_3(ret.a * ret.Te) * ret.sin_wgTe;
  return ret;
}

// B. Doval, C. d'Alessandro and N. Henrich, "The spectrum of glottal flow models",
//   Acta acustica united with acustica, 92(6), 1026-1046, 2006.
FP_TYPE* cig_lfmodel_spectrum(lfmodel model, FP_TYPE* freq, int nf, FP_TYPE* dst_phase) {
  lfparam tmpparam = lfparam_from_lfmodel(model);
  FP_TYPE e = tmpparam.e;
  FP_TYPE a = tmpparam.a;
  FP_TYPE wg = tmpparam.wg;
  FP_TYPE E0 = tmpparam.E0;
  FP_TYPE sin_wgTe = tmpparam.sin_wgTe;
  FP_TYPE cos_wgTe = tmpparam.cos_wgTe;
  FP_TYPE Te = tmpparam.Te;
  FP_TYPE Ta = tmpparam.Ta;
  FP_TYPE T0 = tmpparam.T0;

  FP_TYPE* dst_magn = calloc(nf, sizeof(FP_TYPE));
  for(int i = 0; i < nf; i ++) {
    cplx asubipif = c_cplx(a, - 2.0 * M_PI * freq[i]);
    cplx P1 = c_div(c_cplx(E0, 0), c_add(c_mul(asubipif, asubipif), c_cplx(wg * wg, 0)));
    cplx P2 = c_add(c_cplx(wg, 0),
      c_mul(c_exp(c_mul(asubipif, c_cplx(Te, 0))),
            c_sub(c_mul(asubipif, c_cplx(sin_wgTe, 0)), c_cplx(wg * cos_wgTe, 0))
      ));
    cplx P3 = c_mul(c_cplx(model.Ee, 0),
      c_div(
        c_exp(c_cplx(0, - 2.0 * M_PI * freq[i] * Te)),
        c_mul(c_cplx(0, e * Ta * 2.0 * M_PI * freq[i]), c_cplx(e, 2.0 * M_PI * freq[i]))
      ));
    cplx P4 = c_sub(
      c_mul(c_cplx(e * (1.0 - e * Ta), 0),
            c_sub(c_cplx(1.0, 0),
                  c_exp(c_cplx(0, -2.0 * M_PI * freq[i] * (T0 - Te))))),
            c_cplx(0, e * Ta * 2.0 * M_PI * freq[i]));
    cplx G = c_add(c_mul(P1, P2), c_mul(P3, P4));
    dst_magn[i] = c_abs(G);
    if(isnan(dst_magn[i])) {
      dst_magn[i] = 0;
      G = c_cplx(0, 0);
    }
    if(dst_phase) dst_phase[i] = c_arg(G);
  }
  return dst_magn;
}

FP_TYPE* cig_lfmodel_period(lfmodel model, int fs, int n) {
  lfparam tmpparam = lfparam_from_lfmodel(model);
  FP_TYPE e = tmpparam.e;
  FP_TYPE a = tmpparam.a;
  FP_TYPE wg = tmpparam.wg;
  FP_TYPE E0 = tmpparam.E0;
  FP_TYPE Te = tmpparam.Te;
  FP_TYPE Ta = tmpparam.Ta;
  FP_TYPE T0 = tmpparam.T0;

  int i;
  int ne = round(Te * fs);
  int nT = round(T0 * fs);
  FP_TYPE* y = calloc(n, sizeof(FP_TYPE));
  for(i = 0; i < min(ne, n); i ++) {
    FP_TYPE t = (FP_TYPE)i / fs;
    y[i] = E0 * exp_2(a * t) * sin_2(wg * t);
  }
  for(; i < min(nT, n); i ++) {
    FP_TYPE t = (FP_TYPE)i / fs;
    y[i] = - E0 / e / Ta * (exp_2(- e * (t - Te)) - exp_2(- e * (T0 - Te)));
  }
  return y;
}
