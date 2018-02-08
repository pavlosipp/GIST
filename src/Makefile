########################################################################
#                      Makefile for GIST                               #
########################################################################

########################################################################
#Compiler							       #
########################################################################
#With or without
LIB_NETCDF=no
#using MPI
USE_MPI=no
#Architecture

ifndef HOSTTYPE
  HOSTTYPE := $(shell uname -m)
endif

FC       = ifort #f95i
FFLAGS= -r8 -cpp -O2 -CB -zero -module -I/. $(OBJDIR) 



#using MPI
ifeq ($(USE_MPI),yes)
  FC   = mpif90
  FFLAGS += -DUSE_MPI
endif


ifeq ($(findstring x86_64,${HOSTTYPE}),x86_64)
  LIBS	+= ../lib/libmconf-icc64.a
  FFLAGS += -cxxlib # -pthread

  ifeq ($(LIB_NETCDF),yes)
    LIBS += /afs/ipp/u/stel/netcdf-3.6.2/intel-em64t/lib/libnetcdf.a
    FFLAGS += -DNETCDF
  endif

else
  LIBS += -L/afs/ipp/common/soft/intel/i686/f_9.1/lib
#  LIBS += -lpthread -lcprts -lunwind ../libs/libmconf-icc.a
endif

########################################################################
#Files								       #
########################################################################
include Sources

BINDIR = ../bin
OBJDIR = OBJ
OFILES = $(addprefix $(OBJDIR)/,$(.FFILES:.f=.o))
OFILES1 = $(addprefix $(OBJDIR)/,$(.FFILES:.F90=.o))
MODFILES = $(addprefix $(OBJDIR)/,$(.FFILES:.f=.mod))
MODFILES1 = $(addprefix $(OBJDIR)/,$(.FFILES:.F90=.mod))
DEPENDFILE =   Makedepend

.SUFFIXES:
.SUFFIXES: .f .F90 .o

########################################################################
# Rules								       #
########################################################################

default : install TAGS $(BINDIR)/$(PROJECT)

install :
	@test -d $(BINDIR) || mkdir $(BINDIR)  	
	@test -d $(OBJDIR) || mkdir $(OBJDIR)  	
	@test -d $(INPDIR) || (mkdir $(INPDIR) && cp -r $(TEMPDIR)/input/* $(INPDIR))
#	@test -f ../parameters || cp $(TEMPDIR)/parameters ../.

# Compilation 
$(OBJDIR)/%.o:%.f
	$(FC)  $(FFLAGS) -c -o $@ $<

$(OBJDIR)/%.o:%.F90
	$(FC)  $(FFLAGS) -c -o $@ $<

# Major targets
all : clean install TAGS depend Makefile $(BINDIR)/$(PROJECT) 

clean :
	@rm -f $(OFILES) $(MODFILES,LC) $(BINDIR)/$(PROJECT) $(OBJDIR)/*.mod *.mod #$(DEPENDFILE)
	@rm -f $(OFILES1) $(MODFILES1,LC)
TAGS :
	@etags $(.FFILES)

#Dependencies
depend :
	@makedepf90 -o $(BINDIR)/$(PROJECT) -b $(OBJDIR) $(.FFILES) > $(DEPENDFILE)

$(DEPENDFILE) :
	@makedepf90 -o $(BINDIR)/$(PROJECT) -b $(OBJDIR) $(.FFILES) > $(DEPENDFILE)

include $(DEPENDFILE)
# DO NOT DELETE
