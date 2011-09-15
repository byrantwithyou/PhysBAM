# CUT OFF POINT (lines above will be preserved, lines below will be overwritten)
# Set any necessary USE_* variables (e.g. USE_OPENGL=yes) above CUT OFF POINT.

TARGETS = 3Dmeshing
MAIN_SRC = main.cpp

PUBLIC_LIBRARY_DIR = ../../Public_Library
PERSONAL_LIBRARIES_DIR = ../../Personal_Libraries

LIBRARIES = PhysBAM

include ../../Public_Library/Makefile.common

.PHONY: clean_data
clean_data:
	-rm -f meshing_data/diagnostics/diagnostic* meshing_data/position* meshing_data/velocity* meshing_data/particle* meshing_data/tetrah~
