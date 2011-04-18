# varfrac.mak - Fracture make include file

# Object files #######################################################

# each folder in $VFDIR contains a makefile.include file
# that defines a VF%PATH%OBJ variable containing the list
# of object files it provides

include $(VFDIR)$(S)VF$(S)makefile.include
include $(VFDIR)$(S)Utils$(S)makefile.include

VARFRACOBJ=fisdat$(O) ftdata$(O) fstdout$(O) fwdata$(O) falloc$(O) \
     farray$(O) fplace$(O) ftran$(O) fstep$(O) fivdat$(O) \
     fiadat$(O) fprop$(O) frest$(O) ftransinit$(O) \
     $(VFUTILSOBJ) $(VFOBJ)
		 
#     $(VFOBJDIR)VF_CInterface$(O) $(VFOBJDIR)VF_CInterfacef$(O) \
#     $(VFOBJDIR)VF_Interface$(O) $(VFOBJDIR)fakepmult$(O) $(VFOBJDIR)xdmf$(O)

# Source files #######################################################

#SORC=..$(S)frac$(S)
 
VARFRAC:
	cd $(VFDIR); \
	make FakeVF;     \
	make Utils

