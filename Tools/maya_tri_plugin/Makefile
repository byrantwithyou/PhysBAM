
CXX=/usr/local/compilers/gcc-3.3.2/bin/g++
CXXFLAGS = -g -I. -I.. -I/usr/aw/maya/include -I/usr/X11R6/include -pipe -D_BOOL -DLINUX  -mcpu=pentium4  -Wno-deprecated


all: MAYA_TRI_PLUGIN.so
OBJS=MAYA_TRI_PLUGIN.o MAYA_TRI_PLUGIN_EXTRACTOR.o

MAYA_TRI_PLUGIN.so: ${OBJS}
	${CXX} -shared  -g -pipe -D_BOOL -DLINUX  -mcpu=pentium4  -Wno-deprecated -fno-gnu-keywords -o MAYA_TRI_PLUGIN.so -I. -I.. -I/usr/aw/maya/include -I/usr/X11R6/include ${OBJS} -L/usr/aw/maya/lib -lOpenMaya ../maya_spline_plugin/libPhysBAM_debug.a

