#!/bin/sh
rm -f *.o *.a
gfortran -O3 -c dgedi.f
gfortran -O3 -c dgefa.f
gfortran -O3 -c dpoco.f
gfortran -O3 -c dpofa.f
gfortran -O3 -c dposl.f
gfortran -O3 -c zgedi.f
gfortran -O3 -c zgefa.f

ar rvu liblinpack.a dgedi.o  dgefa.o  dpoco.o  dpofa.o  dposl.o  zgedi.o  zgefa.o
ranlib liblinpack.a
