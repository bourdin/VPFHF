# petsc.mak - Machine and compiler make include file

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

ARCH     = Linux64 
FORT     = ${FC}
LINK     = ${FLINKER} 
CC       = ${CC}
FLEXDIR  = /chap/flexlm/v9.2/libE3-64
CXTC     = $(FLEXDIR)/cxtc.o
#NAN     = ../memman/linux.o
LIBS     = -L $(FLEXDIR) -lpthread

###
### LFLAGS, CFLAGS, FFLAGS are set in petsc makefiles
###

SYSLIB   =
OBJST    = $(GMRSOBJS) $(CXTC)



mpif.h:
	cp /opt/scali/include/mpif.h $(WORK)

$(EXENAM): $(OBJST)
	$(LINK) -o $(EXENAM) $(OBJST)  $(NAN) $(LFLAGS) $(SYSLIB) $(LIBS) $(PETSC_LIB) $(PETSC_FLIB)

clean::
	rm -f $(WORK)/*.f
	rm -f $(WORK)/*.c
	rm -f $(WORK)/*.h
	rm -f $(WORK)/*.mod
	rm -f $(WORK)/*.o
	rm -f $(WORK)/*.i
	rm -f $(WORK)/*.lst
	rm -f $(WORK)/ech
	rm -f $(EXENAM)

