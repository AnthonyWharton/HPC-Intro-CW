CC = cc
COMPILERS := cc gcc icc
OLEVELS := 0 1 2 3

CFLAGS = -std=c99 -Wall
LDFLAGS = -lm
OUTPUT = ./bin/

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

$(OLEVELS): $(OUTPUT)
	$(CC) -O$@ $(CFLAGS) -o $(OUTPUT)jacobi-$(CC)-O$@ jacobi.c $(LDFLAGS)
	@echo $(OUTPUT)jacobi-$(CC)-O$@ >> $(OUTPUT)to-run

_compile_olevels: $(OLEVELS)

$(COMPILERS):
	@$(MAKE) _compile_olevels --no-print-directory CC="$@"

compile-all: $(COMPILERS)

run-all: $(OUTPUT)to-run
	cat $(OUTPUT)to-run | xargs -Iz bash -c 'printf "\n> RUNNING z\n"; z --norder $(JACOBI-NORD) --iterations $(JACOBI-ITER)'

clean:
	rm -rf bin

bluecrystal-job: print-flags compile-all run-all
