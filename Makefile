COMPILER = gcc
GCC = gcc
ICC = icc
COMPILERS := gcc icc
OLEVELS := 2 3

CFLAGS  = -std=c99 -Wall
LDFLAGS = -lm
GCCFLAG = -ffast-math -ftree-vectorizer-verbose=2
ICCFLAG = -xHOST -ipo -no-prec-div -fp-model fast=2 -funroll-loops -qopt-prefetch=4 -mkl=sequential -daal=sequential -qopt-report=5
OUTPUT  = ./bin/

JACOBI-ITER = 20000
JACOBI-NORD = 1000

.NOTPARALLEL: run-all

print-flags:
	@printf "> Flags in operation:\n"
	@printf ">        CC $(CC)\n"
	@printf ">       GCC $(GCC)\n"
	@printf ">       ICC $(ICC)\n"
	@printf ">    CFLAGS $(CFLAGS)\n"
	@printf ">      OPTN $(OPTN)\n"
	@printf ">      OPTH $(OPTH)\n"
	@printf ">   LDFLAGS $(LDFLAGS)\n"
	@printf "\n"

$(OUTPUT):
	mkdir -p $(OUTPUT)

$(OLEVELS):
ifeq ($(COMPILER), $(GCC))
	$(COMPILER) -O$@ $(CFLAGS) $(GCCFLAG) -o $(OUTPUT)jacobi-$(COMPILER)-O$@ jacobi.c $(LDFLAGS)
endif
ifeq ($(COMPILER), $(ICC))
	$(COMPILER) -O$@ $(CFLAGS) $(ICCFLAG) -o $(OUTPUT)jacobi-$(COMPILER)-O$@ jacobi.c $(LDFLAGS)
endif
	@echo $(OUTPUT)jacobi-$(COMPILER)-O$@ >> $(OUTPUT)to-run

_compile_olevels: $(OUTPUT) $(OLEVELS)

$(COMPILERS):
	@$(MAKE) _compile_olevels --no-print-directory COMPILER="$@"

compile-all: $(COMPILERS)

run-all: $(OUTPUT)to-run
	cat $(OUTPUT)to-run | xargs -Iz bash -c 'printf "\n> RUNNING z\n"; z --norder $(JACOBI-NORD) --iterations $(JACOBI-ITER)'

clean:
	rm -rf binr
	rm -rf *.out

bluecrystal-job: print-flags compile-all run-all

bcsubmit:
	bcsubmit -j jacobi.job
