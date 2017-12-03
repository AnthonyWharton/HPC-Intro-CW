################################################################################
###### HPC Coursework 1 Makefile - Anthony Wharton                        ######
################################################################################

### Targets that are safe for use:

# compile-best -  - (DEFAULT) Compiles best config to date
# run-best  -  -  - Compiles and runs the best config to date
# run-all   -  -  - Compiles and runs all configs set in FLAGS
# profile-all  -  - Compiles and profiles all configs set in FLAGS
# compile-all  -  - Compiles all configs set in FLAGS
# compile-all-with-profiling - Same as `compile-all` with profiling support
# clean  -  -  -  - Cleans project directory
# qsub   -  -  -  - Submits the job to BlueCrystal queue
# bluecrystal-job - Target run by the Job File
# bcsubmit  -  -  - Runs bcsubmit - https://github.com/AnthonyWharton/bcsubmit

################################################################################
#################################### FLAGS #####################################
################################################################################
COMPILER = gcc
PROFILING =
GCC = gcc
ICC = icc
COMPILERS := gcc icc
OLEVELS := 2 3

CFLAGS  = -std=c99 -Wall
LDFLAGS = -lm
PFFLAGS = -pg -g
GCCFLAG = -fopenmp
# GCCFLAG = -ffast-math -ftree-vectorizer-verbose=2
ICCFLAG = -qopenmp
# ICCFLAG = -march=native -ipo -no-prec-div -fp-model fast=2 -fp-speculation=fast -funroll-loops -qopt-prefetch=4 -mkl=sequential -daal=sequential -qopt-mem-layout-trans=3 -inline-level=2 -qopt-report=5 -qopenmp
OUTPUT  = ./bin/

JACOBI-ITER = 20000
JACOBI-NORD = 4000

################################################################################
#################################### MISC. #####################################
################################################################################

# Tell make what to *NOT* do in parallel
.NOTPARALLEL: run-all profile-all

# Set Default Target Location
.DEFAULT_GOAL := compile-best

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
######################### GENERIC HIGH LEVEL TARGETS ###########################
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
		gprof -l z >> z.profile'

# Cleans all output files.
clean:
	rm -rvf bin
	rm -rvf *.out
	rm -rvf jacobi

# Bluecrystal Job Target
bluecrystal-job: clean run-all

# Adds job to queue
qsub:
	qsub jacobi.job

# Runs bcsubmit (submits job to bluecrystal from personal machine)
bcsubmit:
	bcsubmit -j jacobi.job

################################################################################
################################### BEST SETUP #################################
################################################################################

# Hardcoded compile flags for best performing Jacobi compilation
compile-best: clean
	icc -O3 -std=c99 -Wall -march=native -ipo -no-prec-div -fp-model fast=2 -fp-speculation=fast -funroll-loops -qopt-prefetch=4 -mkl=sequential -daal=sequential -qopt-mem-layout-trans=3 -inline-level=2 -o jacobi jacobi.c -lm

# Jacobi arguments are still based off the defined variables.
run-best: compile-best
	./jacobi --norder $(JACOBI-NORD) --iterations $(JACOBI-ITER)
