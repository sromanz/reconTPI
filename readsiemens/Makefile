MEX=mex
MKOCTFILE=mkoctfile --mex

CFLAGS=-g -Wall -std=c99
LDFLAGS=-DTEST

all: siemens_meas mexLoadSiemensTraces.oct mexLoadSiemensScandata.oct
	@echo done

clean:
	-rm -r siemens_meas *.o *.dSYM *.mex

%.oct:%.c siemens_meas.o
	CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)" $(MKOCTFILE) -v $^

%.mexa64:%.c
	$(MEX) -v CFLAGS="\$$CFLAGS $(CFLAGS)" LDFLAGS="\$$LDFLAGS $(LDFLAGS)" $^

