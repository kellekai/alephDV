#
# This file will compile the libraries and link to the executable
#

FC = gfortran

FCFLAGS = -std=legacy 
# -fcheck=all -Wall
FLFLAGS = -L${LIBPATH} -lminuit -lquadpack -ltools

LIBPATH = ./lib/
INCPATH = ./inc/
SRCPATH = ./src/

SRC = src/FM.f90 src/FM_CIPT.f90 src/EEM1.f90 src/EEM1_CIPT.f90 \
    src/GEM1.f90 src/GEM1_CIPT.f90 src/FEM.f90 src/FEM_CIPT.f90

OBJ = ${addsuffix .run, ${basename ${notdir ${SRC}}}}

QPSRC = lib/quadpack/dqpsrt.f lib/quadpack/j4save.f lib/quadpack/d1mach.f \
    lib/quadpack/xerhlt.f lib/quadpack/xerprn.f lib/quadpack/fdump.f \
    lib/quadpack/dqagi.f lib/quadpack/xermsg.f lib/quadpack/i1mach.f \
    lib/quadpack/xgetua.f lib/quadpack/dqelg.f lib/quadpack/dqagie.f \
    lib/quadpack/xercnt.f lib/quadpack/dqk15i.f lib/quadpack/gint.f \
    lib/quadpack/xersve.f lib/quadpack/dqag.f lib/quadpack/dqage.f \
    lib/quadpack/dqk21.f lib/quadpack/dqk31.f lib/quadpack/dqk31.f \
    lib/quadpack/dqk51.f lib/quadpack/dqk61.f lib/quadpack/dqk15.f \
    lib/quadpack/dqk41.f lib/quadpack/dgetrf.f lib/quadpack/dgetri.f \
    lib/quadpack/xerbla.f lib/quadpack/ilaenv.f lib/quadpack/dgetf2.f \
    lib/quadpack/dlaswp.f lib/quadpack/dtrsm.f lib/quadpack/dgemm.f \
    lib/quadpack/dtrtri.f lib/quadpack/dgemv.f lib/quadpack/dswap.f \
    lib/quadpack/dlamch.f lib/quadpack/idamax.f lib/quadpack/dscal.f \
    lib/quadpack/dger.f lib/quadpack/lsame.f lib/quadpack/dqk61.f \
    lib/quadpack/dtrmm.f lib/quadpack/dtrti2.f lib/quadpack/dtrmv.f

QPOBJ = ${QPSRC:.f=.o}

MNSRC = lib/minuit/mnrset.F lib/minuit/intract.F lib/minuit/mnparm.F \
    lib/minuit/mnplot.F lib/minuit/mnstat.F lib/minuit/mnpfit.F \
    lib/minuit/mncomd.F lib/minuit/mnderi.F lib/minuit/mnhelp.F \
    lib/minuit/mncuve.F lib/minuit/mnprin.F \
    lib/minuit/mntiny.F lib/minuit/mnsimp.F lib/minuit/mnpsdf.F lib/minuit/mncont.F \
    lib/minuit/mnamin.F lib/minuit/mnunpt.F lib/minuit/mncler.F lib/minuit/mneig.F \
    lib/minuit/mnseek.F lib/minuit/mnread.F lib/minuit/mnsave.F lib/minuit/mnhes1.F \
    lib/minuit/mnexcm.F lib/minuit/mnseti.F lib/minuit/mnmnot.F lib/minuit/mnintr.F \
    lib/minuit/mncros.F lib/minuit/stand.F lib/minuit/mnpout.F lib/minuit/mnbins.F \
    lib/minuit/mninex.F lib/minuit/mninit.F lib/minuit/mnvers.F lib/minuit/mngrad.F \
    lib/minuit/mncrck.F lib/minuit/mnpars.F lib/minuit/mnhess.F lib/minuit/mnrn15.F \
    lib/minuit/mnlims.F lib/minuit/mnset.F lib/minuit/mnmatu.F lib/minuit/mnwerr.F \
    lib/minuit/mnwarn.F lib/minuit/mnpint.F lib/minuit/mnmnos.F lib/minuit/mnline.F \
    lib/minuit/mnfree.F lib/minuit/mnemat.F lib/minuit/mnexin.F lib/minuit/mnerrs.F \
    lib/minuit/mncntr.F lib/minuit/minuit.F lib/minuit/mnfixp.F lib/minuit/mnscan.F \
    lib/minuit/mnrazz.F lib/minuit/mnstin.F lib/minuit/mnmigr.F lib/minuit/mnimpr.F \
    lib/minuit/mninpu.F lib/minuit/mneval.F lib/minuit/mnvert.F lib/minuit/mndxdi.F \
    lib/minuit/mncalf.F

MNOBJ = ${MNSRC:.F=.o}

TLOBJ = complex.o cmplx_root2.o lu.o tipos.o param_dp.o \
    numint.o polint.o cdflib.o

LIBOBJ = lib/libminuit.a lib/libquadpack.a lib/libtools.a 

all: ${OBJ} config clean

FM.run: src/FM.f90 ${LIBOBJ}
	  ${FC} -o $@ $< ${FCFLAGS} ${FLFLAGS}
	
FM_CIPT.run: src/FM_CIPT.f90 ${LIBOBJ}
	  ${FC} -o $@ $< ${FCFLAGS} ${FLFLAGS}
	
EEM1.run: src/EEM1.f90 ${LIBOBJ}
	  ${FC} -o $@ $< ${FCFLAGS} ${FLFLAGS}
	
EEM1_CIPT.run: src/EEM1_CIPT.f90 ${LIBOBJ}
	  ${FC} -o $@ $< ${FCFLAGS} ${FLFLAGS}
	
GEM1.run: src/GEM1.f90 ${LIBOBJ}
	  ${FC} -o $@ $< ${FCFLAGS} ${FLFLAGS}
	
GEM1_CIPT.run: src/GEM1_CIPT.f90 ${LIBOBJ}
	  ${FC} -o $@ $< ${FCFLAGS} ${FLFLAGS}
	
FEM.run: src/FEM.f90 ${LIBOBJ}
	  ${FC} -o $@ $< ${FCFLAGS} ${FLFLAGS}
	
FEM_CIPT.run: src/FEM_CIPT.f90 ${LIBOBJ}
	  ${FC} -o $@ $< ${FCFLAGS} ${FLFLAGS}
	
lib/libminuit.a: ${MNOBJ}
	rm -f $@
	ar rcs $@ $^
	@echo -e "\n***Created Minuit library***\n"

lib/libquadpack.a: ${QPOBJ}
	rm -f $@
	ar rcs $@ $^
	@echo -e "\n***Created Quadpack library***\n"

lib/libtools.a: ${TLOBJ}
	rm -f $@
	ar rcs $@ $^
	@echo -e "\n***Created Tools library***\n"

param_dp.o: lu.o tipos.o
	${FC} -c src/param_dp.f90 
# lu.o tipos.o

numint.o: tipos.o polint.o
	${FC} -c lib/numint.f90 
# tipos.o polint.o

polint.o: tipos.o
	${FC} -c lib/polint.f90 
# tipos.o

lu.o: lib/lu.f90
	${FC} -c lib/lu.f90

tipos.o: lib/tipos.f90
	${FC} -c lib/tipos.f90

cdflib.o: lib/cdflib.f90
	${FC} -c ${FCFLAGS} lib/cdflib.f90

cmplx_root2.o: lib/cmplx_root2.f90
	${FC} -c $<

complex.o: lib/complex.f90 tipos.o param_dp.o
	${FC} -c $<

clean: 
	rm *.mod *.o
	rm lib/quadpack/*.o
	rm lib/minuit/*.o

config:
	cp inc/*.in ./
	
remove: 
	rm ${LIBOBJ} ${OBJ} *.in
