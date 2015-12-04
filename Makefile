######################################
# Makefile written by Zev Kronenberg #
#     zev.kronenberg@gmail.com       #
######################################

CC=gcc
CXX=g++
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always)
CPPFLAGS= -std=c++0x -Wall -DVERSION=\"$(GIT_VERSION)\" 
LIB=-L. -Lvcflib/tabixpp/htslib/
INCLUDE=-I vcflib/include/ -Ivcflib/src/ -Ivcflib/tabixpp/htslib -Ivcflib/tabixpp/
LINKERS=-lhts -lm -lz -lpthread -lsplit

all: mergeSVcallers clean

mergeSVcallers: vcflib/tabixpp/tabix.o libsplit.a libFasta.a libhts.a libvcflib.a
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(LIB)  src/mergeSVcallers.cpp -o mergeSVcallers *.a $(LINKERS) 

libvcflib.a:
	cd vcflib && make && cp libvcflib.a ../
libhts.a: vcflib/libvcflib.a
	cp vcflib/tabixpp/htslib/libhts.a .
libFasta.a: vcflib/libvcflib.a
	ar rcs libfasta.a vcflib/obj/Fasta.o
libsplit.a: vcflib/libvcflib.a
	ar rcs libsplit.a vcflib/src/split.o

clean:
	rm  *.a 