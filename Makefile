CC=gcc
MPICC=mpicc
CILKCC=/usr/local/OpenCilk-9.0.1-Linux/bin/clang
FLAGS=-O1
WARNINGS=-w # tell the compiler to stop emitting warnings.

default: all

sequential:
	$(CC) $(WARNINGS) sequential.c -o sequential head/helpers.c head/mmio.c

pthreads:
	$(CC) $(FLAGS) $(WARNINGS) pthreads.c -o pthreads head/helpers.c head/mmio.c -lpthread

openmp:
	$(MPICC) $(FLAGS) $(WARNINGS) openmp.c -o openmp head/helpers.c head/mmio.c -fopenmp

opencilk:
	$(CILKCC) $(FLAGS) $(WARNINGS) opencilk.c -o opencilk head/helpers.c head/mmio.c -fcilkplus

all: sequential pthreads openmp opencilk

.PHONY: clean

clean:
	rm -f sequential pthreads openmp opencilk

measure_times:
	@printf " ---------- REMAKING DATA.CSV ----------"
	rm -f stats/data.csv
	touch stats/data.csv
	@printf "\n ---------- SEQUENTIAL ----------\n\n"
	./sequential 0
	./sequential 1
	./sequential 2
	./sequential 3
	./sequential 4
	@printf "\n ---------- PTHREAD ----------\n\n"
	./pthreads 0 0 
	./pthreads 0 1
	./pthreads 0 2
	./pthreads 0 3
	./pthreads 0 4
	./pthreads 1 0
	./pthreads 1 1
	./pthreads 1 2
	./pthreads 1 3
	./pthreads 1 4
	./pthreads 2 0
	./pthreads 2 1
	./pthreads 2 2
	./pthreads 2 3
	./pthreads 2 4
	@printf "\n ---------- OPENMP ----------\n\n"
	./openmp 0 0 
	./openmp 0 1
	./openmp 0 2
	./openmp 0 3
	./openmp 0 4
	./openmp 1 0
	./openmp 1 1
	./openmp 1 2
	./openmp 1 3
	./openmp 1 4
	./openmp 2 0
	./openmp 2 1
	./openmp 2 2
	./openmp 2 3
	./openmp 2 4
	@printf "\n ---------- OPENCILK ----------\n\n"
	./opencilk 0 0 
	./opencilk 0 1
	./opencilk 0 2
	./opencilk 0 3
	./opencilk 0 4
	./opencilk 1 0
	./opencilk 1 1
	./opencilk 1 2
	./opencilk 1 3
	./opencilk 1 4
	./opencilk 2 0
	./opencilk 2 1
	./opencilk 2 2
	./opencilk 2 3
	./opencilk 2 4
