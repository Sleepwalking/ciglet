#include "../ciglet.h"

int main(void) {
  FP_TYPE* rand1 = white_noise(1, 10000);
  FP_TYPE* rand2 = white_noise(1, 10000);
  printf("corr = %f, mean = %f, median = %f\n", corr(rand1, rand2, 10000),
    meanfp(rand1, 10000), medianfp(rand1, 10000));
  free(rand1);
  free(rand2);

  int fs, nbit, nx, ny;
  FP_TYPE* x = wavread("test/in2.wav", & fs, & nbit, & nx);
  FP_TYPE* y = rresample(x, nx, 1.5, & ny);
  FP_TYPE* y2 = moving_avg(y, ny, 10);
  FP_TYPE* y2d = diff(y2, ny);
  FP_TYPE* y2c = cumsum(y2d, ny);
  free(y);
  free(y2);
  free(y2d);

  for(int i = 0; i < ny; i ++)
    y2c[i] += randn(0, 0.01 * 0.01);
  wavwrite(y2c, ny, fs * 1.5, nbit, "test/out-resample-mavg.wav");
  free(y2c);

  int nhop = 256;
  int nfft = 2048;
  int nfrm = round(nx / nhop);
  FP_TYPE normfc = 0;
  FP_TYPE** Xm = malloc2d(nfrm, nfft / 2 + 1, sizeof(FP_TYPE));
  FP_TYPE** Xp = malloc2d(nfrm, nfft / 2 + 1, sizeof(FP_TYPE));
  stft(x, nx, nhop, nfrm, 8, 1, & normfc, NULL, Xm, Xp);
  free(x);

  filterbank* mfbank = create_melfilterbank(nfft / 2 + 1, fs / 2, 36, 50, 8000);
  FP_TYPE** Xmfb = filterbank_spgm(mfbank, Xm, nfrm, nfft, fs, 0);
  FP_TYPE** Xmfcc = calloc(nfrm, sizeof(FP_TYPE*));
  for(int i = 0; i < nfrm; i ++)
    Xmfcc[i] = be2cc(Xmfb[i], 36, 12, 0);

  figure* fg = plotopen();
  imagesc(fg, Xmfcc, nfrm, 12);
  plotclose(fg);
  FP_TYPE** Xmmf = calloc(nfft / 2 + 1, sizeof(FP_TYPE*));
  FP_TYPE** Xmtr = transpose(Xm, nfrm, nfft / 2 + 1, sizeof(FP_TYPE));
  for(int i = 0; i < nfft / 2 + 1; i ++) {
    Xmmf[i] = medfilt1(Xmtr[i], nfrm, 25);
    for(int j = 0; j < nfrm; j ++)
      Xmmf[i][j] = log_2(Xmmf[i][j]);
  }
  free2d(Xmtr, nfft / 2 + 1);
  Xmtr = transpose(Xmmf, nfft / 2 + 1, nfrm, sizeof(FP_TYPE));
  free2d(Xmmf, nfft / 2 + 1);
  for(int i = 0; i < nfrm; i ++) {
    FP_TYPE* tmp = medfilt1(Xmtr[i], nfft / 2 + 1, 25);
    free(Xmtr[i]);
    Xmtr[i] = tmp;
  }
  free2d(Xmtr, nfrm);

  free2d(Xmfb, nfrm);
  free2d(Xmfcc, nfrm);
  delete_filterbank(mfbank);

  FP_TYPE** Xm2 = copy2d(Xm, nfrm, nfft / 2 + 1, sizeof(FP_TYPE));
  for(int i = 0; i < nfrm; i ++) {
    for(int j = round(30 + pow(sin(i * 0.05), 2) * 200); j < nfft / 2 + 1; j ++)
      Xm2[i][j] = 0;
  }

  y = istft(Xm2, Xp, nhop, nfrm, 8, 1, normfc, & ny);
  wavwrite(y, ny, fs, nbit, "test/out-stft-wow.wav");
  free(y);
  free2d(Xm2, nfrm);

  FP_TYPE** C = malloc2d(nfrm, nfft / 2 + 1, sizeof(FP_TYPE));
  for(int i = 0; i < nfrm; i ++) {
    free(spec2env(Xm[i], nfft, fs, 250.0, C[i]));
  }
  FP_TYPE** Xm3 = cegm2spgm(C, nfrm, nfft, nfft / 2 + 1);
  fg = plotopen();
  imagesc(fg, Xm3, nfrm, nfft / 2 + 1);
  plotclose(fg);

  FP_TYPE* Xm3_flattened = flatten(Xm3, nfrm, nfft / 2 + 1, sizeof(FP_TYPE));
  FP_TYPE Xm3median = medianfp(Xm3_flattened, nfrm * (nfft / 2 + 1));
  FP_TYPE Xm3mean = meanfp(Xm3_flattened, nfrm * (nfft / 2 + 1));
  printf("Median of spectrogram: %f\n", Xm3median);
  printf("Mean of spectrogram: %f (stdvar: %f)\n", Xm3mean,
    sqrt(varfp(Xm3_flattened, nfrm * (nfft / 2 + 1))));
  //for(int i = 0; i < nfrm * (nfft / 2 + 1); i ++)
  //  Xm3_flattened[i] = exp_2(Xm3_flattened[i]);
  FP_TYPE** Xm4 = reshape(Xm3_flattened, nfrm, nfft / 2 + 1, sizeof(FP_TYPE));
  free(Xm3_flattened);
  free2d(Xm3, nfrm); free2d(Xm4, nfrm); free2d(C, nfrm);

  free2d(Xm, nfrm); free2d(Xp, nfrm);
  return 0;
}
