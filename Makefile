################################################################################
###### HPC Coursework 1 Makefile - Anthony Wharton                        ######
################################################################################

################################################################################
#################################### FLAGS #####################################
################################################################################
COMPILER = gcc
PROFILING = ""
GCC = gcc
ICC = icc
COMPILERS := icc
OLEVELS := 3

CFLAGS  = -std=c99 -Wall
LDFLAGS = -lm
PFFLAGS = -pg -g
GCCFLAG = -ffast-math -ftree-vectorizer-verbose=2
ICCFLAG = -xHOST -ipo -no-prec-div -fp-model fast=2 -fp-speculation=fast -funroll-loops -qopt-prefetch=4 -mkl=sequential -daal=sequential -qopt-mem-layout-trans=3 -inline-level=2 -qopt-report=5
OUTPUT  = ./bin/

JACOBI-ITER = 20000
JACOBI-NORD = 2000

################################################################################
############################ MISC. TARGETS #####################################
################################################################################

# Tell make what to *NOT* do in parallel
.NOTPARALLEL: run-all profile-all

# Creates the output location
$(OUTPUT):
	mkdir -p $(OUTPUT)

################################################################################
############################# ACTUAL COMPILATION ###############################
################################################################################

# Compiles the program with the compiler specified in $(COMPILER) from the
# recursive make calls in the $(COMPILERS) target, OR with any profiling
# support added by compile-all-with-profiling into $(PROFILING) with  all the
# O-Levels specified in $(OLEVELS) (referenced with $@ in the target), into the
# $(OUTPUT) location.
$(OLEVELS):
ifeq ($(COMPILER), $(GCC))
	$(COMPILER) -O$@ $(CFLAGS) $(PROFILING) $(GCCFLAG) -o $(OUTPUT)jacobi-$(COMPILER)-O$@ jacobi.c $(LDFLAGS)
else
ifeq ($(COMPILER), $(ICC))
	$(COMPILER) -O$@ $(CFLAGS) $(PROFILING) $(ICCFLAG) -o $(OUTPUT)jacobi-$(COMPILER)-O$@ jacobi.c $(LDFLAGS)
else
  $(COMPILER) -O$@ $(CFLAGS) $(PROFILING)            -o $(OUTPUT)jacobi-$(COMPILER)-O$@ jacobi.c $(LDFLAGS)
endif
endif
	@echo $(OUTPUT)jacobi-$(COMPILER)-O$@ >> $(OUTPUT)to-run

################################################################################
############################# SET UP COMPILATION ###############################
################################################################################

# Hidden target used by $(COMPILERS)
_compile_olevels: $(OUTPUT) $(OLEVELS)
# Recursively calls make to compile the program with the different compilers
# specifying the compilers in $(COMPILERS), referenced by $@ below.
$(COMPILERS):
	@$(MAKE) _compile_olevels --no-print-directory COMPILER="$@"

# Compiles all programs with all compilers using the $(COMPILERS) target.
compile-all: $(COMPILERS)

# Hidden target used by the compile-all-with-profiling target
_compile_profiling: $(COMPILERS)
# Recursively calls make with specified $(PROFILING) flags specified
# in $(PFFLAGS)
compile-all-with-profiling:
	@$(MAKE) _compile_profiling  --no-print-directory PROFILING="$(PFFLAGS)"

################################################################################
############################# HIGH LEVEL TARGETS ###############################
################################################################################

# Compiles the program with all compilers specified in $(COMPILERS) and with
# O-Levels specified in $(OLEVELS), and then runs all programs.
run-all: compile-all $(OUTPUT)to-run
	@cat $(OUTPUT)to-run | xargs -Iz bash -c ' \
		printf "\n> RUNNING z\n"; \
		z --norder $(JACOBI-NORD) --iterations $(JACOBI-ITER);'

# Compiles the program with all compilers specified in $(COMPILERS) and with
# O Levels specified in $(OLEVELS), and then runs all programs through gprof.
profile-all: compile-all-with-profiling $(OUTPUT)to-run
	cat $(OUTPUT)to-run | xargs -Iz bash -c ' \
		printf "\n> PROFILING z\n"; \
		z --norder $(JACOBI-NORD) --iterations $(JACOBI-ITER); \
		gprof -l z >> $(OUTPUT)/profile'

# Cleans all output files.
clean:
	rm -rvf bin
	rm -rvf *.out

# Bluecrystal Job Target
bluecrystal-job: compile-all run-all

# Runs bcsubmit (submits job to bluecrystal from personal machine)
bcsubmit:
	bcsubmit -j jacobi.job
