# GNUmakefile 
# Created for MaGe modules by Jason Detwiler, June 2005

include ../buildTools/config.mk

name := $(notdir $(shell pwd))
G4TARGET := $(name)
G4EXLIB := true

.PHONY: all
all: #lib  -> nothing to build in analysis for now...

include /home/admin/geant4.10.04/config/binmake.gmk

LIBFILES := \
	libROOT_ioROOT.so libROOT_ioROOT.rootmap libROOT_ioROOT_rdict.pcm \
	libROOT_mjioROOT.so libROOT_mjioROOT.rootmap libROOT_mjioROOT_rdict.pcm

install:
	mkdir -p $(DESTDIR)$(libdir)
	@for f in $(LIBFILES); do \
	  echo $(INSTALL_DATA) $$f $(DESTDIR)$(libdir)/; \
	  $(INSTALL_DATA) $$f $(DESTDIR)$(libdir)/; \
	done

uninstall:
	@for f in $(LIBFILES); do \
	  echo $(RM) $(DESTDIR)$(libdir)/`basename $$f`; \
	  $(RM) $(DESTDIR)$(libdir)/`basename $$f`; \
	done
