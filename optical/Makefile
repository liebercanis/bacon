include $(MGDODIR)/buildTools/config.mk

# Give the list of applications, which must be the stems of cc files with 'main'.  There
# can be more than one.  In our example, this means there is a test.cc and a test1.cc
#APPS = anaOptical anaBACoN anaMapBACoN anaMap anaOptical anaRadio binHistogram mapMerger readNtuple anaCosmic
APPS =  anaOptical  anaRadio anaCosmic 
# The next two lines are important
LOCAL=/home/gold/XenonDoping
SHLIB =
ARCHIVE =
TAMDIR = /home/admin/MGDO/tam
# Include the correct flags,
INCLUDEFLAGS = $(CLHEP_INCLUDE_FLAGS) -I$(MGDODIR)/Base -I$(MGDODIR)/Root -I$(MGDODIR)/Transforms -I$(LOCAL)/io -I$(G4INSTALL)/include/Geant4/

INCLUDEFLAGS += -I$(MGDODIR)/Majorana -I$(MGDODIR)/MJDB $(ROOT_INCLUDE_FLAGS) -I$(TAMDIR)/inc -I$(TAMDIR)/include -I$(MGDODIR)/Tabree

LIBFLAGS = -L$(MGDODIR)/lib -lMGDORoot -lMGDOBase -lMGDOTransforms
LIBFLAGS += -L/usr/local/lib -lCLHEP-Random-2.4.1.0
#-lMGDOMajorana -lMGDOGerdaTransforms -lMGDOMJDB -lMGDOTabree
LIBFLAGS += $(ROOT_LIB_FLAGS) -lSpectrum -L$(TAMDIR)/lib -lTAM
LIBFLAGS += -L/usr/local/root/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTDataFrame -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic

include $(MGDODIR)/buildTools/BasicMakefile

.PHONEY: show cleanup

show:
	@echo TAMDIR $(TAMDIR)
	@echo MGDODIR $(MGDODIR)
	@echo " "
	@echo $(INCLUDEFLAGS)
	@echo " "
	@echo $(LIBFLAGS)

cleanup:
	@echo "CLEANUP"
	rm *.o anaOptical anaBACoN  anaMap anaRadio anaCosmic
