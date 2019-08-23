PREFIX = /usr
CC = $(CROSS)gcc
LINK = $(CROSS)gcc
AR = $(CROSS)ar

ifeq 'Darwin' '$(shell uname)'
  CFLAGS_PLAT =
else
  CFLAGS_PLAT = -fopenmp
endif

CFLAGS = -Iexternal -D_POSIX_C_SOURCE=2 -DFP_TYPE=float -Og -g -std=c99 -Wall -fPIC $(CFLAGSEXT) $(CFLAGS_PLAT)
ARFLAGS = -rv
OBJS = ciglet.o fftsg.o fastmedian.o wavfile.o

default: libciglet.a

test: ciglet-test
	./ciglet-test noplot

test-fft: ciglet-test-fft
	./ciglet-test-fft

single-file: single-file/ciglet.c single-file/ciglet.h
	$(CC) $(CFLAGS) -o single-file/ciglet.o -c single-file/ciglet.c -Wno-unused-result

ciglet-test: libciglet.a test/test.c
	$(CC) $(CFLAGS) test/test.c libciglet.a -o ciglet-test -lm

ciglet-test-fft: libciglet.a test/test-fft.c
	$(CC) $(CFLAGS) test/test-fft.c libciglet.a -o ciglet-test-fft -lm

libciglet.a: $(OBJS)
	$(AR) $(ARFLAGS) libciglet.a $(OBJS)
	@echo Done.

install: libciglet.a ciglet.h
	mkdir -p $(PREFIX)/lib/ $(PREFIX)/include/
	cp libciglet.a $(PREFIX)/lib
	cp ciglet.h $(PREFIX)/include
	cp external/fastapprox-all.h $(PREFIX)/include

ciglet.o: ciglet.c ciglet.h
fftsg.o: external/fftsg_h.c
	$(CC) $(CFLAGS) -o fftsg.o -c external/fftsg_h.c
fastmedian.o: external/fast_median.c
	$(CC) $(CFLAGS) -o fastmedian.o -c external/fast_median.c
wavfile.o: external/wavfile.c
	$(CC) $(CFLAGS) -o wavfile.o -c external/wavfile.c -Wno-unused-result

%.o: %.c
	$(CC) $(CFLAGS) -o $*.o -c $*.c

single-file/ciglet.c: ciglet.c
	mkdir -p single-file
	cat external/fftsg_h.c > single-file/ciglet.c
	cat external/fast_median.c >> single-file/ciglet.c
	cat external/wavfile.c >> single-file/ciglet.c
	cat ciglet.c >> single-file/ciglet.c

single-file/ciglet.h: ciglet.h
	mkdir -p single-file
	cat external/fastapprox-all.h > single-file/ciglet.h
	echo "" >> single-file/ciglet.h
	echo "#define CIGLET_SINGLE_FILE" >> single-file/ciglet.h
	echo "" >> single-file/ciglet.h
	cat ciglet.h >> single-file/ciglet.h

clean:
	@echo 'Removing all temporary binaries... '
	@rm -f *.o *.a
	@echo Done.

clear:
	@echo 'Removing all temporary binaries... '
	@rm -f *.o *.a
	@echo Done.
