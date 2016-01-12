#
# --- HYCOM 2.0 makefile 
#
# --- Tunable parameters in ../config/$(ARCH)_$(TYPE)
#

.SUFFIXES: 
.SUFFIXES: .c .F .f .o

.F:
	@echo "Must have an explicit rule for" $*
.f:
	@echo "Must have an explicit rule for" $*
.c:
	@echo "Must have an explicit rule for" $*

include ../config/$(ARCH)_$(TYPE)

default: hycom


MODS =   mod_xc.o mod_za.o mod_pipe.o

OBJS =	 hycom.o archiv.o  barotp.o bigrid.o  blkdat.o cnuity.o convec.o \
	diapfl.o dpthuv.o  dpudpv.o forfun.o  geopar.o hybgen.o \
	icloan.o inicon.o  inikpp.o  inimy.o  latbdy.o \
	matinv.o momtum.o  mxkprf.o  mxkrt.o  mxkrtm.o  mxpwp.o \
	overtn.o poflat.o  prtmsk.o  psmoo.o restart.o stencl.o \
	thermf.o trcupd.o  tsadvc.o  zebra.o \
       machine.o  wtime.o machi_c.o

SRCS =	 hycom.f archiv.f  barotp.f bigrid.f  blkdat.f cnuity.f convec.f \
	diapfl.f dpthuv.f  dpudpv.f forfun.f  geopar.f hybgen.f \
	icloan.f inicon.f  inikpp.f  inimy.f  latbdy.f \
	matinv.f momtum.f  mxkprf.f  mxkrt.f  mxkrtm.f  mxpwp.f \
	overtn.f poflat.f  prtmsk.f  psmoo.f restart.f stencl.f \
	thermf.f trcupd.f  tsadvc.f  zebra.f \
       machine.f  wtime.f machi_c.c \
      mod_pipe.F \
        mod_xc.F mod_xc_sm.F mod_xc_mp.F \
        mod_za.F mod_za_sm.F mod_za_mp.F mod_za_mp1.F

hycom:	$(MODS) $(OBJS)
	$(LD)  $(LDFLAGS) -o hycom  $(MODS) $(OBJS) $(EXTRALIBS)

hycom.o:   hycom.f   mod_xc.o common_blocks.h stmt_fns.h mod_za.o mod_pipe.o
archiv.o:  archiv.f  mod_xc.o common_blocks.h            mod_za.o
barotp.o:  barotp.f  mod_xc.o common_blocks.h                     mod_pipe.o
bigrid.o:  bigrid.f  mod_xc.o 
blkdat.o:  blkdat.f  mod_xc.o common_blocks.h
cnuity.o:  cnuity.f  mod_xc.o common_blocks.h                     mod_pipe.o
convec.o:  convec.f  mod_xc.o common_blocks.h stmt_fns.h
diapfl.o:  diapfl.f  mod_xc.o common_blocks.h stmt_fns.h
dpthuv.o:  dpthuv.f  mod_xc.o common_blocks.h
dpudpv.o:  dpudpv.f  mod_xc.o 
forfun.o:  forfun.f  mod_xc.o common_blocks.h            mod_za.o
geopar.o:  geopar.f  mod_xc.o common_blocks.h stmt_fns.h mod_za.o
hybgen.o:  hybgen.f  mod_xc.o common_blocks.h stmt_fns.h          mod_pipe.o
icloan.o:  icloan.f  mod_xc.o common_blocks.h stmt_fns.h
inicon.o:  inicon.f  mod_xc.o common_blocks.h stmt_fns.h          mod_pipe.o
inikpp.o:  inikpp.f  mod_xc.o common_blocks.h stmt_fns.h
inimy.o:   inimy.f   mod_xc.o common_blocks.h stmt_fns.h
latbdy.o:  latbdy.f  mod_xc.o common_blocks.h
machine.o: machine.F
machi_c.o: machi_c.c
matinv.o:  matinv.f  mod_xc.o common_blocks.h
momtum.o:  momtum.f  mod_xc.o common_blocks.h stmt_fns.h          mod_pipe.o
mxkprf.o:  mxkprf.f  mod_xc.o common_blocks.h stmt_fns.h
mxkrt.o:   mxkrt.f   mod_xc.o common_blocks.h stmt_fns.h
mxkrtm.o:  mxkrtm.f  mod_xc.o common_blocks.h stmt_fns.h
mxpwp.o:   mxpwp.f   mod_xc.o common_blocks.h stmt_fns.h
overtn.o:  overtn.f  mod_xc.o common_blocks.h
poflat.o:  poflat.f
prtmsk.o:  prtmsk.f
psmoo.o:   psmoo.f   mod_xc.o 
restart.o: restart.f mod_xc.o common_blocks.h            mod_za.o
stencl.o:  stencl.f  mod_xc.o common_blocks.h
thermf.o:  thermf.f  mod_xc.o common_blocks.h stmt_fns.h
trcupd.o:  trcupd.f  mod_xc.o common_blocks.h
tsadvc.o:  tsadvc.f  mod_xc.o common_blocks.h stmt_fns.h          mod_pipe.o
wtime.o:   wtime.F
zebra.o:   zebra.f   mod_xc.o 

mod_pipe.o: \
         mod_pipe.F  mod_xc.o common_blocks.h
mod_xc.o:  mod_xc.F  dimensions.h mod_xc_sm.F mod_xc_mp.F
mod_za.o:  mod_za.F  mod_xc.o     mod_za_sm.F mod_za_mp.F mod_za_mp1.F

tar: Makefile README.* *.h $(SRCS)
	tar -cvf hycom_src.tar Makefile README.* *.h $(SRCS)
