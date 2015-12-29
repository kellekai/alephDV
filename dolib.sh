#!/bin/bash
gfortran -c i1mach.f
gfortran -c xgetua.f
gfortran -c fdump.f
gfortran -c j4save.f
gfortran -c xercnt.f
gfortran -c xerprn.f
gfortran -c xersve.f
gfortran -c xerhlt.f
gfortran -c xermsg.f
gfortran -c d1mach.f
gfortran -c dqpsrt.f
gfortran -c dqk15i.f
gfortran -c dqelg.f
gfortran -c dqagie.f
gfortran -c dqagi.f
gfortran -c gint.f
gfortran -c tipos.f90
gfortran -c lu.f90
gfortran -c polint.f90 tipos.o
gfortran -c numint.f90 tipos.o polint.o
gfortran -c param_dp.f90 lu.o tipos.o
ar rcs libchi2fov.a *.o *.mod
mv ./libchi2fov.a ./lib/libchi2fov.a
rm *.o


