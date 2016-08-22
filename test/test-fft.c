#include "../ciglet.h"
#include <sys/time.h>

static double gettime() {
  struct timeval t;
  gettimeofday(&t, NULL);
  return (t.tv_sec + (t.tv_usec / 1000000.0)) * 1000.0;
}

int main(void) {
  int N = 512;
  FP_TYPE* in_re  = calloc(N, sizeof(FP_TYPE));
  FP_TYPE* in_im  = calloc(N, sizeof(FP_TYPE));
  FP_TYPE* out_re = calloc(N, sizeof(FP_TYPE));
  FP_TYPE* out_im = calloc(N, sizeof(FP_TYPE));
  FP_TYPE* buffer = calloc(N * 2, sizeof(FP_TYPE));

  double T = 0;
  for(int i = 0; i < 1000; i ++) {
    for(int j = 0; j < N; j ++) {
      #if 0
        in_re[j] = randu();
        in_im[j] = randu();
      #else
        in_re[j] = 1e-44; // subnormal for single precision floats
        in_im[j] = -1e-44;
      #endif
    }
    double t0 = gettime();
    ifft(in_re, in_im, out_re, NULL, N, buffer);
    double t1 = gettime();
    T += t1 - t0;
  }
  printf("Time taken: %fms\n", T);

  free(in_re); free(in_im); free(out_re); free(out_im);
}
