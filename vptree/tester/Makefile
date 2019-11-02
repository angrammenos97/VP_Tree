# define the shell to bash
SHELL := /bin/bash

# define the C/C++ compiler to use,default here is clang
CC = gcc-7

all: test_sequential test_pthreads test_openmp test_cilk

test_sequential:
	tar -xvzf code.tar.gz
	cd vptree; make lib; cd ..
	cd vptree; cp lib/*.a inc/vptree.h ../; cd ..
	$(CC) tester.c vptree_sequential.a -o $@ -lm
	./test_sequential

test_pthreads:
	tar -xvzf code.tar.gz
	cd vptree; make lib; cd ..
	cd vptree; cp lib/*.a inc/vptree.h ../; cd ..
	$(CC) -pthread tester.c vptree_pthreads.a -o $@ -lm
	./test_pthreads

test_openmp:
	tar -xvzf code.tar.gz
	cd vptree; make lib; cd ..
	cd vptree; cp lib/*.a inc/vptree.h ../; cd ..
	$(CC) -omp tester.c vptree_openmp.a -o $@ -lm
	./test_openmp

test_cilk:
	tar -xvzf code.tar.gz
	cd vptree; make lib; cd ..
	cd vptree; cp lib/*.a inc/vptree.h ../; cd ..
	$(CC) -fcilkplus tester.c vptree_cilk.a -o $@ -lm
	./test_cilk

clean:
	$(RM) -fr ./vptree vptree* test_*
