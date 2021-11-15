CC=gcc
MPICC=mpicc
CILKCC=/usr/local/OpenCilk-9.0.1-Linux/bin/clang
FLAGS=-O1
WARNINGS=-w # tell the compiler to stop emitting warnings.

default: all

sequential:
	$(CC) $(WARNINGS) sequential.c -o sequential head/helpers.c head/mmio.c

pthreads:
	$(CC) $(WARNINGS) pthreads.c -o pthreads head/helpers.c head/mmio.c -lpthread

openmp:
	$(MPICC) $(FLAGS) $(WARNINGS) openmp.c -o openmp head/helpers.c head/mmio.c -fopenmp

opencilk:
	$(CILKCC) $(FLAGS) $(WARNINGS) opencilk.c -o opencilk head/helpers.c head/mmio.c -fcilkplus

all: sequential pthreads openmp opencilk

.PHONY: clean

clean:
	rm -f sequential pthreads openmp opencilk

test:
	@printf "\n ---------- SEQUENTIAL ----------\n\n"
	./sequential
	@printf "\n ---------- PTHREAD ----------\n\n"
	./pthreads
	@printf "\n ---------- OPENMP ----------\n\n"
	./openmp
	@printf "\n ---------- OPENCILK ----------\n\n"
	./opencilk
