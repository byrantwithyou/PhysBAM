# Run "make depend" (and later, if you suspect or know file dependencies have changed, if make gives you weird errors...)
# Run "make" to create your project in release mode, and "make debug" to make a debugging executable
# Run "make clean" to clean up (e.g. before changing from debugging to performance flags)

################################ edit the following for your project and system

# where the PhysBAM files reside (edit this if needed)
PHYSBAM_LIB=../../Public_Library

# the name of the executable (edit this)
MAIN_PROGRAM = subdivide_tri_file

# .cpp source files for your project (edit this)
MAIN_SRC = main.cpp
ADVECTION_EQUATIONS_SRC = ADVECTION_HAMILTON_JACOBI_WENO.cpp ADVECTION.cpp
ARRAYS_SRC =
BOUNDARIES_SRC = BOUNDARY_CONSTANT_EXTRAPOLATION.cpp
DEFORMABLE_OBJECTS_SRC = 
FORCES_AND_TORQUES_SRC = 
GEOMETRY_SRC = TRIANGLE_3D.cpp TRIANGULATED_SURFACE.cpp PLANE.cpp BOX_3D.cpp TRIANGLE_HIERARCHY.cpp SEGMENT_3D.cpp SPHERE.cpp OCTREE_VOLUME.cpp
GRIDS_SRC = TRIANGLE_MESH.cpp SEGMENT_MESH.cpp TRIANGLE_SUBDIVISION.cpp OCTREE_MESH.cpp
LEVEL_SETS_SRC = OCTREE_LEVELSET.cpp
MATH_TOOLS_SRC =
MATRICES_AND_VECTORS_SRC =
NONLINEAR_EQUATIONS_SRC =
OPENGL_SRC =
PARTICLES_SRC = PARTICLE_3D.cpp PARTICLE.cpp
RANDOM_NUMBERS_SRC =
READ_WRITE_SRC =
RENDERING_SRC =
RIGID_BODIES_SRC =

# compiler specifics (edit this)
#
# for example:
#   CC = g++ -Wall -march=i686
#   RELEASE_FLAGS = -O3 -fomit-frame-pointer -DNDEBUG
#   DEBUG_FLAGS = -g 
#   LINK = g++
#   RELEASE_LINKFLAGS =
#   DEBUG_LINKFLAGS =
#   DEPEND_OUTPUT_FLAG = -MT    # this is for gcc 3 and up
#   #DEPEND_OUTPUT_FLAG = -o    # this is for obsolete gcc 2.9x (e.g. cygwin, default RedHat install)
#
include Makefile.defs

# here you can append extra flags to RELEASE_FLAGS, LINK_FLAGS, etc.
# but please make sure they are platform-independent!
LINK += -lm

################################## the rest of the file shouldn't need editing

include ../../Public_Library/Makefile.common

